from pandas import DataFrame
import distance_util as du
import pandas as pd
import numpy as np
from BaseClasses import ResidueFeatureExtractor,ChainFeatureExtractor,SeqFeatureExtractor
from Bio.PDB.Entity import Entity
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Data import *

from typing import Iterable, List,Tuple,Union,Dict,Any

from Bio import BiopythonParserWarning
from util import _integrated_residue_iterator,allowed_residue_source,allowd_scheme,_impute_default_value
from distance_util import residue_distance_matrix
from CdrAnnotator import _run_anarci,_impute_cdr,_impute_vicinity
from BiopyInternalFeature import _impute_sasa
import warnings




def _impute_static_feature_struct(object:allowed_residue_source):
    for residue in _integrated_residue_iterator(object):
        resname:str=residue.get_resname()
        residue.xtra['hydrophobicity']=hydrophobicity_scale.get(resname,0)
        residue.xtra['static_charge']=charge_scale.get(resname,0)
        residue.xtra['i_sasa']=SASA_scale.get(resname,np.mean(list(SASA_scale.values())))
        residue.xtra['georged_i_sasa']=GeorgeDSASA_scale.get(resname,np.mean(list(SASA_scale.values())))

def _impute_salt_bridge(object:allowed_residue_source):
    if len(object.child_list)>1:
        warnings.warn("multi-frame object detected."
             "only frame 0 will be processed.",
             BiopythonParserWarning,)
    positive_name=['NZ','NH1','NH2']
    negative_name=['OD1','OD2','OE1','OE2']
    positive_atoms=[i for i in object.get_atoms() if i.get_name() in positive_name]
    negative_atoms=[i for i in object.get_atoms() if i.get_name() in negative_name]
    for positive_atom in positive_atoms:
        vicinity_dict=du.atom_within_threshold(negative_atoms,positive_atom,3.2)
        if len(vicinity_dict)>0:
            positive_atom.get_parent().xtra['salt_bridge']=1
            for negative_atom in vicinity_dict.keys():
                negative_atom.get_parent().xtra['salt_bridge']=1
    for residue in _integrated_residue_iterator(object):
        # if residue.xtra.get('salt_bridge',0)==1:
        #     print(residue.get_full_id(),end='\t')
        # else:
        #     residue.xtra['salt_bridge']=0
        if residue.xtra.get('salt_bridge',0)!=1:
            residue.xtra['salt_bridge']=0

def _impute_salted_charge_n_hydro(object:allowed_residue_source):
    '''
    must impute static_feature and salt_bridge first
    '''
    for residue in _integrated_residue_iterator(object):
        residue.xtra['salted_charge'] = residue.xtra['static_charge'] if not residue.xtra['salt_bridge'] else 0
        residue.xtra['salted_hydrophobicity'] = residue.xtra['hydrophobicity'] if not residue.xtra['salt_bridge'] else hydrophobicity_scale['GLY']

def _impute_fsasa(object:allowed_residue_source):
    '''
    must run after `_impute_static_feature_struct` and `_impute_sasa` 
    '''
    for residue in _integrated_residue_iterator(object):
        residue.xtra['f_sasa']=0.01*residue.xtra['SASA_INTERNAL']/residue.xtra['i_sasa']

def _impute_psh_like_feature(distance_matrix:pd.Series,input_key:str,output_key:str,func=lambda x:x,threshhold:float=7.5)->None:
    '''
    distance_matrix:must comes from du.residue_distance_matrix
    threshhold<=0 means no threshhold
    '''
    for row in distance_matrix.iteritems():
        # assert (input_key in row[0][0].xtra and 
        #                 input_key in row[0][1].xtra),f'input key not ready in {str(row[0][0])} or {str(row[0][1])}'
        if row[1]==0:
            continue
        else:
            if threshhold>0:
                value=func(row[0][0].xtra[input_key])*func(row[0][1].xtra[input_key])/row[1]**2 if row[1]<threshhold else 0
            else:
                value=func(row[0][0].xtra[input_key])*func(row[0][1].xtra[input_key])/row[1]**2
        row[0][0].xtra[output_key]=row[0][0].xtra.get(output_key,0)+value/4
        row[0][1].xtra[output_key]=row[0][1].xtra.get(output_key,0)+value/4

def _impute_psh(distance_matrix:pd.Series)->None:
    _impute_psh_like_feature(distance_matrix,input_key='salted_hydrophobicity',output_key='psh',func=lambda x:(x+13.5)/9,threshhold=7.5)

def _impute_ppc(distance_matrix:pd.Series)->None:
    _impute_psh_like_feature(distance_matrix,input_key='salted_charge',output_key='ppc',func=lambda x:x if x>0 else 0 ,threshhold=7.5)

def _impute_pnc(distance_matrix:pd.Series)->None:
    _impute_psh_like_feature(distance_matrix,input_key='salted_charge',output_key='pnc',func=lambda x:x if x<0 else 0 ,threshhold=7.5)


def _impute_surface(object:allowed_residue_source,fsasa_threshold=0.075):
    '''
    must run `_impute_fsasa` first.
    '''
    for residue in _integrated_residue_iterator(object):
        if residue.xtra['f_sasa']>0.075:
            residue.xtra['surface']=1
        else:
            residue.xtra['surface']=0

def _impute_surface_charge(object:allowed_residue_source):
    '''
    must run `_impute_fsasa` and `_impute_salted_charge` first.
    '''
    for residue in _integrated_residue_iterator(object):
        residue.xtra['surface_charge']=residue.xtra['salted_charge']*residue.xtra['surface']

class StaticPropertyImputer(ResidueFeatureExtractor):
    def __init__(self)->None:
        ResidueFeatureExtractor.__init__(self,operation_name='static',canonical_only=True)
    
    def _produce_feature(self):
        
        pass


class static_feature_imputer(SeqFeatureExtractor):
    def _produce_feature(self):
        hydrophobicity_list=[]
        static_charge_list=[]
        i_sasa_list=[]
        georged_i_sasa_list=[]
        for res in self.seq:
            resname=amino1to3dict[res]
            hydrophobicity_list.append(hydrophobicity_scale[resname])
            static_charge_list.append(charge_scale[resname])
            i_sasa_list.append(SASA_scale[resname])
            georged_i_sasa_list.append(GeorgeDSASA_scale[resname])
        self.frame['hydrophobicity']=hydrophobicity_list
        self.frame['static_charge']=static_charge_list
        self.frame['i_sasa']=i_sasa_list
        self.frame['georged_i_sasa']=georged_i_sasa_list


class TAPExtractor(ResidueFeatureExtractor):
    def __init__(self,scheme:allowd_scheme='a'):
        self.scheme=scheme
        ResidueFeatureExtractor.__init__(self,operation_name='extendedTAP',canonical_only=True)

    def _produce_feature(self):
        #impute antibody specific property:
        #Fv_fragment_id(Fv_annotation) CDR_id(CDR) vicinity chain_type(H,L,K)
        self._produce_sequence()
        _run_anarci(self.sequences,self.scheme,self.object)
        _impute_cdr(self.object,self.scheme)
        self.vicinity_residue=_impute_vicinity(self.object)

        #impute basic residue property:
        #SASA_INTERNAL,i_sasa,f_sasa,surface
        #static_charge,salt_bridge,salted_charge
        #

        _impute_sasa(self.object)
        _impute_static_feature_struct(self.object)
        _impute_salt_bridge(self.object)
        _impute_salted_charge_n_hydro(self.object)
        _impute_fsasa(self.object)
        _impute_surface(self.object)
        _impute_surface_charge(self.object)

        #impute psh ppc pnc
        self.vicinity_distance_matrix=residue_distance_matrix(self.vicinity_residue)
        _impute_psh(self.vicinity_distance_matrix)
        _impute_default_value(self.object,'psh',0)
        _impute_ppc(self.vicinity_distance_matrix)
        _impute_default_value(self.object,'ppc',0)
        _impute_pnc(self.vicinity_distance_matrix)
        _impute_default_value(self.object,'pnc',0)

        self._object_feature_to_frame()
        self._reduce()

    def _reduce(self):
        '''
        self.frame(residue-level feature pool to Fv-level feature pool)
        '''
        fv_fragments=self.frame[self.frame['Fv_annotation'].apply(lambda x: True if 'Fv' in x else False)].groupby(['chain','Fv_annotation'])
        self.fv_frame=pd.DataFrame()
        
        def _count_cdr(iter:Iterable)->int:
                sum=0
                for i in iter:
                    if 'CDR'in i:
                        sum +=1
                return sum
        
        self.fv_frame['chain_type']=fv_fragments['chain_type'].apply(lambda x :x.iloc[0])
        self.fv_frame['CDR_lenth']=fv_fragments['CDR'].apply(_count_cdr)
        self.fv_frame['PSH']=fv_fragments['psh'].sum()
        self.fv_frame['PPC']=fv_fragments['ppc'].sum()
        self.fv_frame['PNC']=fv_fragments['pnc'].sum()
        self.fv_frame['sc_SFvCSP']=fv_fragments['surface_charge'].sum()

class ExtendedTAPExtractor(ResidueFeatureExtractor):
    def __init__(self)->None:
        ResidueFeatureExtractor.__init__(self,operation_name='extendedTAP',canonical_only=True)
    def _produce_feature(self):

        if not self._CDR_vicinity:
            pass
