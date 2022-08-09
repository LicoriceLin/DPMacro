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
from .util import integrated_residue_iterator,write_out
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

def model(residues:allowed_residue_source)->Model:
    residue_list=list(integrated_residue_iterator(residues))
    residue_list.sort(key=lambda x :-ord(x.get_parent().id)*10000+x.id[1])
    model=Model('0')
    for residue in residue_list:
        chainid=residue.get_parent().id
        if chainid not in model:
            new_chain=Chain(chainid)
            model.add(new_chain)
        model[chainid].add(residue.copy())
    return model

class PocketExtractor(ResidueFeatureExtractor):
# class structure_aligner:
    def __init__(self):
        '''
        set rna=True if processing RNA
        so far it only affects the generation of self.sequences and self.reference_sequences
        '''
        ResidueFeatureExtractor.__init__(self,operation_name='AFill_aligner', 
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
        self._find_proteins()

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
                        if i.id[0]==' ']
    
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
        self.pocket=residue_within_threshold(self.proteins,self.selected_ligand,distance)
        # self.neighbor_chain= tuple(set([i.get_parent().id for i in self.neighbor]))
        # add return

    def set_chain_as_ligand(self,chain:Chain):
        assert chain in list(self.reference.get_chains()),'invalid chain'
        self.selected_ligand=chain
        _other_residue=[ residue for residue in self.proteins if residue not in list(chain.get_residues()) ]
        self.neighbor=residue_within_threshold(_other_residue,self.selected_ligand,6)

    def save(self,prefix:str=""):
        '''
        save the ligand and the pocket with the given profix 
        '''
        write_out(model(self.selected_ligand),f'{prefix}_ligand.pdb')
        write_out(model(self.pocket),f'{prefix}_pocket.pdb')

    def transform(self,**kwargs) :
        raise NotImplementedError
       
class PocketFilter:
    pass
        
