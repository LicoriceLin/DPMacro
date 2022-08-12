from typing import Union,Literal,Dict,Tuple,List
import tempfile
from itertools import permutations
import pandas as pd

from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

import pymol

from .BaseClasses import ResidueFeatureExtractor
from .util import Structure,Residue,Chain,Model
from .util import RESIDUE_HOLDER,allowed_residue_source
from .util import integrated_residue_iterator,write_out,model
from .distance_util import residue_within_threshold 
from .Data import amino3to1dict,nucleic_acid_dict

   

allowed_platform=Literal['pymol','vmd']
def reslist_to_guicode(reslist:allowed_residue_source,guiplatform:allowed_platform='pymol')->str:
    '''
    generate the selection sentence for pymol / vmd.
    note: 
    it only generate the sub-model level selection.
    that means,for example, you load xxx.pdb into pymol and want to select part of it by this function
    the selection should be f"xxx and { reslist_to_guicode(...) }" 
    '''
    def chain_n_residue(residue:Residue)->str:
        chain=residue.parent
        if guiplatform =='pymol':
            return '(chain ' + str(chain.id)  + ' and resi ' + str(residue.id[1])+' )'
        elif guiplatform=='vmd':
            return  '(chain ' + str(chain.id) + ' and resid ' + str(residue.id[1])+' )'
    return '( '+' or '.join([chain_n_residue(i) for i in integrated_residue_iterator(reslist) 
                                            if i.parent.id != 'place_holder'])+' )'

class PocketExtractor(ResidueFeatureExtractor):
# class structure_aligner:
    def __init__(self):
        '''
        set rna=True if processing RNA
        so far it only affects the generation of self.sequences and self.reference_sequences
        '''
        ResidueFeatureExtractor.__init__(self,operation_name='PocketExtractor', 
                                            canonical_only=False)

    def load_structure(self,struct:Union[str,Structure])->None:
        '''
        public method of self._set_structure()
        1.it read in a structure as origin `_set_structure` does.
        2.it will init self.sequence,a {chainid:seq} dict, according to structure.
            i.X(unrecognised residue e.g. solvent) at the begining and end of the sequence will be removed
            ii.protein-rna hybrid system is untested. if process RNA, remember to set rna=True when init the instance.
        3.it will init self.object_path,a path for pymol to readin input structure.
        '''
        self._set_structure(struct)
        assert len(self.object)==1,'please plit the frames first!'
        self._find_ligands()
        # self._find_proteins()

    def _find_ligands(self):
        '''
        set the self.ligands, a list (may transfer to a Tuple later) containing all ligand Residue instance in obj.
        solvent will be removed by built-in rules of biopython.
        '''
        self.ligands:List[Residue]=[i for i in integrated_residue_iterator(self.object) 
                        if i.id[0] not in [' ','W']]

    def _find_proteins(self):
        '''
        set the self.protein,a list (may transfer to a Tuple later) containing all protein Residue instance in obj.
        '''
        self.proteins:List[Residue]=[i for i in integrated_residue_iterator(self.object) 
                        if i not in self.ligands and i.id[0]!='W']
    
    def set_ligand(self,ligand:Residue,distance:float=6):
        '''
        choose a ligand from self.ligands
        set it to self.selected_ligand
        calculate protein residues within 6 A of self.selected_ligand
        set the self.neighbor as a list of these residues.
        set the self.neighbor_chain as a list of chain contains self.neighbor.
        '''
        assert ligand in self.ligands,'invalid ligand!'
        self.selected_ligand=ligand
        self.pocket=residue_within_threshold([residue for residue in integrated_residue_iterator(self.object) if not residue is self.selected_ligand ],self.selected_ligand,distance)
        self._to_model()
        # self.neighbor_chain= tuple(set([i.get_parent().id for i in self.neighbor]))
        # add return

    def set_chain_as_ligand(self,chain:Chain,distance:float=6):
        assert chain in list(self.reference.get_chains()),'invalid chain'
        self.selected_ligand=chain
        _other_residue=[ residue for residue in self.proteins if residue not in list(chain.get_residues()) ]
        self.pocket=residue_within_threshold(_other_residue,self.selected_ligand,distance)
        self._to_model()

    def flexible_set_ligand(self,reses:allowed_residue_source,distance:float=6):
        self.selected_ligand=reses
        _other_residue=[ residue for residue in self.object.get_residues() if residue not in list(integrated_residue_iterator(reses)) ]
        self.pocket=residue_within_threshold(_other_residue,self.selected_ligand,distance)
        self._to_model()
    def _to_model(self):
        self.ligand_model=model(self.selected_ligand)
        self.pocket_model=model(self.pocket)
        self.lig_poc_model=model(list(integrated_residue_iterator(self.selected_ligand))+self.pocket)
    def save(self,prefix:str=""):
        '''
        save the ligand and the pocket with the given profix 
        '''
        write_out(self.ligand_model,f'{prefix}_ligand.pdb')
        write_out(self.pocket_model,f'{prefix}_pocket.pdb')

    def save_as_whole(self,out='ligand-complex.pdb'):
        
        write_out(self.lig_poc_model,out)
    def transform(self,**kwargs) :
        raise NotImplementedError
       
class PocketFilter:
    pass
        
