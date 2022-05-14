
import pandas as pd
import tempfile
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from typing import Union,Literal,Set

from BaseClasses import ResidueFeatureExtractor
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from util import RESIDUE_HOLDER,allowed_residue_source
from util import read_in,_integrated_residue_iterator,sequence_from_object
from distance_util import atom_within_threshold 

   

allowed_platform=Literal['pymol','vmd']
def reslist_to_guicode(reslist:allowed_residue_source,guiplatform:allowed_platform)->str:
    def chain_n_residue(residue:Residue)->str:
        chain=residue.parent
        if guiplatform =='pymol':
            return '(chain ' + str(chain.id)  + ' and resi ' + str(residue.id[1])+' )'
        elif guiplatform=='vmd':
            return  '(chain ' + str(chain.id) + ' and resid ' + str(residue.id[1])+' )'
    return ' or '.join([chain_n_residue(i) for i in _integrated_residue_iterator(reslist) 
                                            if i.parent.id != 'place_holder'])

def write_fasta(prefix:str,seq_dict:dict,chainid:Union[str,None]=None,output:str='poc.fasta'):
    with open(output,'a') as f:
        if chainid is None:
            for chain,seq in seq_dict.items():
                f.write('>'+prefix+chain+'\n')
                f.write(seq+'\n')
        else:
            f.write('>'+prefix+chainid+'\n')
            f.write(seq_dict[chainid]+'\n')

def residue_within_threshold(entity1:allowed_residue_source,entity2:allowed_residue_source,threshold)->Set[Residue]:
    '''
    '''
    return set(_integrated_residue_iterator(
                            atom_within_threshold(entity1,entity2,threshold).keys()))



class structure_aligner(ResidueFeatureExtractor):
# class structure_aligner:
    def __init__(self):
        ResidueFeatureExtractor.__init__(self,operation_name='AFill_aligner', 
                                            canonical_only=False)
    def _set_structure(self,struct:Union[str,Structure])->None:
        '''
        this method is reload!
        '''
        if isinstance(struct,str):
            self.object=read_in(struct,struct)
        elif isinstance(struct,Structure):
            self.object=struct.copy()
        else:
            raise TypeError
        self.sequences=sequence_from_object(self.object)
    def _set_reference_structure(self,struct:Union[str,Structure])->None:
        '''
        '''
        if isinstance(struct,str):
            self.reference=read_in(struct,struct)
        elif isinstance(struct,Structure):
            self.reference=struct.copy()
        else:
            raise TypeError
        self.reference_sequences=sequence_from_object(self.reference)
    def _find_ligands(self):
        self.ligands=[i for i in _integrated_residue_iterator(self.reference) 
                        if i.id[0] not in [' ','W']]
    def _find_proteins(self):
        self.proteins=[i for i in _integrated_residue_iterator(self.reference) 
                        if i.id[0]==' ']
    def _align(self,object_chain:str,reference_chain:str)->pd.DataFrame:
        with tempfile.TemporaryDirectory() as dirname:
            fasta_file=dirname+'/align.fasta'
            aln_file=dirname+'/align.aln'
            write_fasta(prefix=self.object.id+'_',seq_dict=self.sequences,chainid=object_chain,output=fasta_file)
            write_fasta(prefix=self.reference.id+'_',seq_dict=self.reference_sequences,chainid=reference_chain,output=fasta_file)
            clustalw_cline=ClustalwCommandline("clustalw2", infile=fasta_file,outfile=aln_file)
            clustalw_cline()
            self.alignment = AlignIO.read(aln_file,format='clustal')
        pointer=0
        iter1=self.object[0][object_chain].get_residues()
        iter2=self.reference[0][reference_chain].get_residues()
        res_list1=[]
        res_list2=[]
        while pointer<len(self.alignment._records[0].seq):
            record1=self.alignment._records[0].seq[pointer]
            record2=self.alignment._records[1].seq[pointer]
            if record1!='-':
                residue1=next(iter1)
            else:
                residue1=RESIDUE_HOLDER
            if record2!='-':
                residue2=next(iter2)
            else:
                residue2=RESIDUE_HOLDER

            res_list1.append(residue1)
            res_list2.append(residue2)

            pointer +=1

        self.alignment_frame=pd.DataFrame()
        self.alignment_frame['object_residue']=res_list1
        self.alignment_frame['reference_residue']=res_list2

    def map_neighbor(self):
        f=self.alignment_frame
        if not hasattr(self,'mapped_neighbor'):
            self.mapped_neighbor=f[f['reference_residue'].isin(self.neighbor)]['object_residue'].to_list()
        else:
            self.mapped_neighbor.extend(f[f['reference_residue'].isin(self.neighbor)]['object_residue'].to_list())
    def load_reference(self,struct:Union[str,Structure]):
        self._set_reference_structure(struct)
        self._find_ligands()
        self._find_proteins()

    def set_ligand(self,ligand:Residue):
        assert ligand in self.ligands,'invalid ligand!'
        self.selected_ligand=ligand
        self.neighbor=residue_within_threshold(self.proteins,self.selected_ligand,6)

    def _produce_feature(self,object_chain:str,reference_chain:str):
        self._align(object_chain,reference_chain)

    def transform(self, struct: Union[str, Structure]) -> pd.DataFrame:
        '''
        '''
        raise NotImplementedError
        
        
