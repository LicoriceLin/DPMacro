from typing import Dict
import tempfile
import os

from .BaseClasses import ResidueFeatureExtractor
from Bio.PDB.Structure import Structure
from Bio.PDB import DSSP,ShrakeRupley,HSExposureCB

from .util import write_out
from .util import integrated_residue_iterator,allowed_residue_source,impute_default_value

def impute_dssp(object:Structure,dssp_argsargs:Dict={}):
    '''
    run dssp on input structure.
    allow all extra args to initiate DSSP instance.
    just input them in a dict as the second parameter.

    impute these feature in input's Residue.xtra dict:
        secondary structure information:
        'SS_DSSP':str

        conformation:
        'PHI_DSSP': float(degree),
        'PSI_DSSP': float(degree),
        
        
        
        energy and hydrogen bond information:
        (Maybe I need to check the detail for each of them LOL)
        'NH_O_1_RELIDX_DSSP'
        'NH_O_1_ENERGY_DSSP'
        'O_NH_1_RELIDX_DSSP'
        'O_NH_1_ENERGY_DSSP'
        'NH_O_2_RELIDX_DSSP'
        'NH_O_2_ENERGY_DSSP'
        'O_NH_2_RELIDX_DSSP'
        'O_NH_2_ENERGY_DSSP'
        
        SASA:
        'EXP_DSSP_RASA': float,
        'EXP_DSSP_ASA':float

        (an invalid 'feature':
        'DSSP_INDEX': int,residue index )

    Note:
    only full structure instance is allowed as input
    only protein-only system is tested
    
    '''
    with tempfile.TemporaryDirectory() as dirname:
        _pdb=os.path.join(dirname,'tmp.pdb')
        write_out(object,_pdb)
        DSSP(object[0],_pdb,**dssp_argsargs)
    impute_default_value(object,'SS_DSSP','x')
    for key in ['NH_O_1_RELIDX_DSSP',
        'NH_O_1_ENERGY_DSSP',
        'O_NH_1_RELIDX_DSSP',
        'O_NH_1_ENERGY_DSSP',
        'NH_O_2_RELIDX_DSSP',
        'NH_O_2_ENERGY_DSSP',
        'O_NH_2_RELIDX_DSSP',
        'O_NH_2_ENERGY_DSSP',
        'EXP_DSSP_RASA',
        'EXP_DSSP_ASA',
        'PHI_DSSP',
        'PSI_DSSP',
        'DSSP_INDEX']:
        
        impute_default_value(object,key,-1)


def impute_sasa(object:Structure,sasa_argsargs:Dict={}):
    '''
    run biopython's SASA module on the input structure.
    allow all extra args to initiate ShrakeRupley instance.
    just input them in a dict as the second parameter.
    
    impute 'SASA_INTERNAL' in input's Residue.xtra dict.
        warning:this value can vary quit a lot from SASA in dssp.
    
    Note:
    allow all Entity, but only recommend Structure and Model. Other input may be unstable.
    '''
    _ShrakeRupley=ShrakeRupley(**sasa_argsargs)
    _ShrakeRupley.compute(object,level="R")
    for residue in object.get_residues():
        residue.xtra['SASA_INTERNAL']=residue.sasa
        del residue.sasa

def impute_hse(object:Structure,hse_argsargs:Dict={}):
    '''
    calculate the half-sphere exposure value use biopython's internal function, 
    which is a new index to show if a residue is buried.

    allow all extra args to initiate HSExposureCB instance.
    just input them in a dict as the second parameter.

    impute these feature in input's Residue.xtra dict:
    'EXP_HSE_B_U': int,number of CA upper along the CA-CB axis
    'EXP_HSE_B_D': int,number of CA downside along the CA-CB axis
    'EXP_HSE_RATIO': float,ratio of the 2 value below

    Note:
    only full structure instance is allowed as input
    only protein-only system is tested

    '''
    HSExposureCB(object,**hse_argsargs)
    impute_default_value(object,'EXP_HSE_B_U',-1)
    impute_default_value(object,'EXP_HSE_B_D',1)
    for residue in object.get_residues():
        residue.xtra['EXP_HSE_RATIO']=residue.xtra['EXP_HSE_B_U']/residue.xtra['EXP_HSE_B_D'] if residue.xtra['EXP_HSE_B_D'] != 0 else 20

def impute_angles(object:allowed_residue_source):
    '''
    calculate angle & dihedral in the input structure.

    impute these feature in input's Residue.xtra dict:
    'PHI_INTERNAL','PSI_INTERNAL','OMG_INTERNAL','TAU_INTERNAL','CHI1_INTERNAL'

    chi2 or more is ok for impute. but not implement right now

    Note:
    all input form in util.allowed_residue_source is allowed.
    but ligands and non-canonical component in protein is not supported right now.

    '''
    object.atom_to_internal_coordinates()
    for residue in integrated_residue_iterator(object):
        _ric=residue.internal_coord.get_angle
        ### the default value may be better for other value than 0?
        residue.xtra['PHI_INTERNAL']=_ric('phi') if _ric('phi') else 0
        residue.xtra['PSI_INTERNAL']=_ric('psi') if _ric('psi') else 0
        residue.xtra['OMG_INTERNAL']=_ric('omg') if _ric('omg') else 0
        residue.xtra['TAU_INTERNAL']=_ric('tau') if _ric('tau') else 0
        residue.xtra['CHI1_INTERNAL']=_ric('chi1') if _ric('chi1') else 0


class BiopyInternalFeature(ResidueFeatureExtractor):
    '''
    subclass of `ResidueFeatureExtractor`.
    a very explicit wrapper to run all above at a time for a structure.
    initiate it with args you need, and `transform` pdb file or structures one after another. 
    (e.g. if you'd like to run internal sasa; what's the args for its execution.)
    '''
    def __init__(self,dssps:bool=True,sasa:bool=False,hse:bool=True,angles:bool=True,
                    dssp_argsargs:Dict={},sasa_args:Dict={},hse_args:Dict={}):
        '''
        '''
        self._dssps=dssps
        self._dssp_args=dssp_argsargs
        self._sasa=sasa
        self._sasa_args=sasa_args
        self._hse=hse
        self._hse_args=hse_args
        self._angles=angles
        ResidueFeatureExtractor.__init__(self,operation_name='BiopyInternal',canonical_only=True)

    def _produce_feature(self):
        '''
        '''
        if self._dssps:
            # with tempfile.TemporaryDirectory() as dirname:
            #     _pdb=os.path.join(dirname,'tmp.pdb')
            #     write_out(self.object,_pdb)
            #     DSSP(self.object[0],_pdb,**self._dssp_args)
            impute_dssp(self.object,self._dssp_args)
        
        if self._sasa:
            # _ShrakeRupley=ShrakeRupley(**self._sasa_args)
            # _ShrakeRupley.compute(self.object,level="R")
            # for _residue in self.object.get_residues():
            #     _residue.xtra['SASA_INTERNAL']=_residue.sasa
            #     del _residue.sasa
            impute_sasa(self.object,self._sasa_args)

        if self._hse:
            impute_hse(self.object,self._hse_args)
            # HSExposureCB(self.object[0],**self._hse_args)
            # for _residue in self.object.get_residues():
            #     _residue.xtra['EXP_HSE_RATIO']=_residue.xtra['EXP_HSE_B_U']/_residue.xtra['EXP_HSE_B_D']
            
        if self._angles:
            impute_angles(self.object)
            # self.object.atom_to_internal_coordinates()
            # for _residue in self.object.get_residues():
            #     _ric=_residue.internal_coord.get_angle
            #     ### the default value may be better for other value than 0?
            #     _residue.xtra['PHI_INTERNAL']=_ric('phi') if _ric('phi') else 0
            #     _residue.xtra['PSI_INTERNAL']=_ric('psi') if _ric('psi') else 0
            #     _residue.xtra['OMG_INTERNAL']=_ric('omg') if _ric('omg') else 0
            #     _residue.xtra['TAU_INTERNAL']=_ric('tau') if _ric('tau') else 0
            #     _residue.xtra['CHI1_INTERNAL']=_ric('chi1') if _ric('chi1') else 0

        self._object_feature_to_frame()


if __name__=='__main__':
    import sys
    testpdb=sys.argv[1]
    Bpi=BiopyInternalFeature(dssp_argsargs={'acc_array':"Wilke"},
                            sasa_args={'probe_radius':1.40, 'n_points':300},
                            hse_args={'radius':18, 'offset':0})
    Bpi.transform(testpdb)
    if not os.path.isdir('test'):
        os.mkdir('test')
    Bpi.frame.to_csv('test/testBiopyInternalFeature.csv')
