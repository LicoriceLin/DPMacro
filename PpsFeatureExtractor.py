from typing import Iterable,Callable,Dict
import copy
import tempfile

from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from numpy import matrix
from BaseClasses import ResidueFeatureExtractor
from distance_util import residue_distance_matrix,residue_within_threshold
from util import allowed_residue_source,integrated_residue_iterator,sequence_from_object,impute_default_value
from Data import hydrophobic_set,Charg_set,Aroma_set,SASA_scale,amino1to3dict
from BiopyInternalFeature import impute_sasa,impute_angles,impute_hse,impute_dssp
from TAPFeature import impute_surface,impute_static_feature,impute_fsasa
from TAPFeature import impute_psh,impute_ppc,impute_pnc,impute_salt_bridge,impute_salted_charge_n_hydro

def impute_PXX(object:allowed_residue_source,key:str,ref_set:Iterable[str]):
    '''
    a trivial feature in PremPS:
    if a residue belongs to a certain residue-type set (charge, aroma, etc.) 
    and is on the preotein surface, PXX will be set to 1; else 0.
    
    you can define surface by TAPFeature.impute_surface, 
    or define a function on you own, which impute 1/0 into Residue.xtra['surface']

    must run `impute_surface` first
    in PremPS, threshold of surface should be set to ( fSASA > )0.2 
    in TAP, it's 0.075.

    such function is used to describe the protein system (calculate how many xx residue on the surface) in PPS, 
    rather than the mutation site.
    '''
    for residue in integrated_residue_iterator(object):
        if residue.resname in ref_set and residue.xtra['surface']:
            residue.xtra[key]=1
        else:
            residue.xtra[key]=0

def impute_PFWY(object:allowed_residue_source):
    '''
    see `impute_PXX`
    '''
    impute_PXX(object,'PFWY',Aroma_set)

def impute_PRKDE(object:allowed_residue_source):
    '''
    see `impute_PXX`
    '''
    impute_PXX(object,'PRKDE',Charg_set)

def impute_PL(object:allowed_residue_source):
    '''
    see `impute_PXX`
    '''
    impute_PXX(object,'PL',{'LEU',})

def impute_PLIV(object:allowed_residue_source):
    '''
    see `impute_PXX`
    '''
    impute_PXX(object,'PLIV',{'LEU','ILE','VAL'})

def impute_Nxx(object:Structure,key:str,ref_set:Iterable[str],neighbor_range:int=12):
    '''
    a trivial feature in PremPS:
    count how many residues of a certain type (charge, aroma, etc.) 
    is within `neighbor_range` of each residue.

    define residue to count in `ref_set`,
    define the key to impute in Residue.xtra in `key`
    define the range of neighbor in `neighbor_range` ( e.g. 12 -> left 11 + object residue + right 11)

    such function is used to describe the mutation site in PPS.

    '''
    assert len(list(object.get_models()))==1,'please split multi-frame file by util.split_multiframe first.'
    sequences=sequence_from_object(object)
    for chainid,chain_seq in sequences.items():
        # print(len(chain_seq))
        range_list=list(range(len(chain_seq)))
        residue_list=list(object[0][chainid].get_residues())
        for res, loc_id in zip(residue_list,range_list):
            neighbor_seq=chain_seq[max(0,loc_id-neighbor_range+1):min(len(chain_seq),loc_id+neighbor_range)]
            # print(neighbor_seq,end='\t')
            three_code_neighbor_seq=[ amino1to3dict[residue]for residue in neighbor_seq]
            Nxx=sum([three_code_neighbor_seq.count(i) for i in ref_set])
            res.xtra[key]=Nxx
            # print(Nxx)

def impute_NHydro(object:Structure):
    '''
    see `impute_Nxx`
    '''
    impute_Nxx(object,'NHydro',hydrophobic_set)
       
def impute_NCharg(object:Structure):
    '''
    see `impute_Nxx`
    '''
    impute_Nxx(object,'NCharg',Charg_set)

def impute_pssm(object:Structure):
    assert len(list(object.get_models()))==1,'please split multi-frame file by util.split_multiframe first.'
    sequences=sequence_from_object(object)
    for chainid,chain_seq in sequences.items():
        with tempfile.TemporaryDirectory() as dirname:
            pass

#used 
def impute_sNxx(object:Structure,key:str,ref_set:Iterable[str],neighbor_range:int=6):
    '''
    structure version of Nxx,search for spacial vicinity rather than sequence vicinity.   
    '''
    for residue in integrated_residue_iterator(object):
        neighbor_list=[residue.resname for residue in residue_within_threshold(object,residue)]
        Nxx=sum([neighbor_list.count(i) for i in ref_set])
        residue.xtra[key]=Nxx



class PPSExtractor(ResidueFeatureExtractor):
    def __init__(self):
        ResidueFeatureExtractor.__init__(self,operation_name='PPS',canonical_only=True)
    def _produce_feature(self):

        #impute_internal_feature
        impute_dssp(self.object)
        impute_angles(self.object)
        impute_hse(self.object)
        impute_sasa(self.object)

        #impute surface
        impute_static_feature(self.object)
        impute_salt_bridge(self.object)
        impute_salted_charge_n_hydro(self.object)
        impute_fsasa(self.object)
        impute_surface(self.object,fsasa_threshold=0.2)

        #impute_Pxx
        impute_PFWY(self.object)
        impute_PRKDE(self.object)
        impute_PL(self.object)
        impute_PLIV(self.object)

        #impute_Nxx
        impute_NHydro(self.object)
        impute_NCharg(self.object)

        #tentative:impute extened_pxx_of_tap
        impute_salt_bridge(self.object)
        impute_salted_charge_n_hydro(self.object)
        self.surface_distance_matrix=residue_distance_matrix([i for i in self.object.get_residues() if i.xtra['surface']==1 ])
        impute_psh(self.surface_distance_matrix)
        impute_ppc(self.surface_distance_matrix)
        impute_pnc(self.surface_distance_matrix)
        [impute_default_value(self.object,i,0) for i in ['ppc','pnc','psh']]

        self._object_feature_to_frame()

    def reduce(self,residue:Residue):
        assert residue in list(self.object.get_residues()),'invalid residue'
        output_dict={}
        for key in [#basic description
                    'hydrophobicity','static_charge',
                    
                    'EXP_DSSP_ASA','EXP_DSSP_RASA','SASA_INTERNAL','i_sasa','georged_i_sasa','f_sasa',
                    'surface','salt_bridge','salted_charge','salted_hydrophobicity',
                    
                    #DSSP description
                    'SS_DSSP','NH_O_1_RELIDX_DSSP','NH_O_1_ENERGY_DSSP','O_NH_1_RELIDX_DSSP',
                    'O_NH_1_ENERGY_DSSP','NH_O_2_RELIDX_DSSP','NH_O_2_ENERGY_DSSP',
                    'O_NH_2_RELIDX_DSSP','O_NH_2_ENERGY_DSSP',

                    #angles
                    'PHI_INTERNAL','PSI_INTERNAL','OMG_INTERNAL','TAU_INTERNAL','CHI1_INTERNAL',
                    
                    #hse
                    'EXP_HSE_B_U','EXP_HSE_B_D','EXP_HSE_RATIO',
                    
                    #pps
                    'NHydro','NCharg',
                    
                    #TAP-like feature
                    'psh','ppc','pnc']:
            output_dict[key]=residue.xtra[key]
        
        chain_belong=residue.get_parent().id
        for key in ['PFWY','PRKDE','PL','PLIV']:
            output_dict[key]=self.frame.groupby('chain')[key].sum()[chain_belong]
        return output_dict


