from BaseClasses import ResidueFeatureExtractor
from Bio.PDB.Structure import Structure
import typing
import tempfile
import os
from util import write_out
from Bio.PDB import DSSP,ShrakeRupley,HSExposureCB
class BiopyInternalFeature(ResidueFeatureExtractor):
    '''
    '''
    def __init__(self,dssps:bool=True,sasa:bool=False,hse:bool=True,angles:bool=True,
                    dssp_argsargs:typing.Dict={},sasa_args:typing.Dict={},hse_args:typing.Dict={}):
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
            with tempfile.TemporaryDirectory() as dirname:
                _pdb=os.path.join(dirname,'tmp.pdb')
                write_out(self.object,_pdb)
                DSSP(self.object[0],_pdb,**self._dssp_args)
        
        if self._sasa:
            _ShrakeRupley=ShrakeRupley(**self._sasa_args)
            _ShrakeRupley.compute(self.object,level="R")
            for _residue in self.object.get_residues():
                _residue.xtra['SASA_INTERNAL']=_residue.sasa
                del _residue.sasa

        if self._hse:
            HSExposureCB(self.object[0],**self._hse_args)
            for _residue in self.object.get_residues():
                _residue.xtra['EXP_HSE_RATIO']=_residue.xtra['EXP_HSE_B_U']/_residue.xtra['EXP_HSE_B_D']
            
        if self._angles:
            self.object.atom_to_internal_coordinates()
            for _residue in self.object.get_residues():
                _ric=_residue.internal_coord.get_angle
                ### the default value may be better for other value than 0?
                _residue.xtra['PHI_INTERNAL']=_ric('phi') if _ric('phi') else 0
                _residue.xtra['PSI_INTERNAL']=_ric('psi') if _ric('psi') else 0
                _residue.xtra['OMG_INTERNAL']=_ric('omg') if _ric('omg') else 0
                _residue.xtra['TAU_INTERNAL']=_ric('tau') if _ric('tau') else 0
                _residue.xtra['CHI1_INTERNAL']=_ric('chi1') if _ric('chi1') else 0

        self._object_feature_to_frame()