from typing import Union,Literal,Dict,Tuple,List
import tempfile
from itertools import permutations
import pandas as pd

from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

import pymol

from .BaseClasses import ResidueFeatureExtractor
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from .util import RESIDUE_HOLDER,allowed_residue_source
from .util import read_in,integrated_residue_iterator,sequence_from_object,write_fasta,wash_sequence
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

class structure_aligner(ResidueFeatureExtractor):
# class structure_aligner:
    def __init__(self,rna:bool=False):
        '''
        set rna=True if processing RNA
        so far it only affects the generation of self.sequences and self.reference_sequences
        '''
        ResidueFeatureExtractor.__init__(self,operation_name='AFill_aligner', 
                                            canonical_only=False)
        self.resname_dict=nucleic_acid_dict if rna else  amino3to1dict
    def _set_structure(self,struct:Union[str,Structure],filepath:Union[str,None]=None)->None:
        '''
        this method is reload.
        1.it read in a structure as origin `_set_structure` does.
        2.it will init self.sequence,a {chainid:seq} dict, according to structure.
            i.X(unrecognised residue e.g. solvent) at the begining and end of the sequence will be removed
            ii.protein-rna hybrid system is untested. if process RNA, remember to set rna=True when init the instance.
        3.it will init self.object_path,a path for pymol to readin input structure.
        '''
        if isinstance(struct,str):
            self.object=read_in(struct,struct)
            self.object_path=struct
        elif isinstance(struct,Structure):
            assert isinstance(filepath,str),'must give me the path of of pdb'
            self.object=struct.copy()
            self.object_path=filepath
        else:
            raise TypeError
        dirty_seqs=sequence_from_object(self.object,self.resname_dict)
        self.sequences=wash_sequence(dirty_seqs)

    def load_structure(self,struct:Union[str,Structure],filepath:Union[str,None]=None)->None:
        '''
        public method of self._set_structure()
        1.it read in a structure as origin `_set_structure` does.
        2.it will init self.sequence,a {chainid:seq} dict, according to structure.
            i.X(unrecognised residue e.g. solvent) at the begining and end of the sequence will be removed
            ii.protein-rna hybrid system is untested. if process RNA, remember to set rna=True when init the instance.
        3.it will init self.object_path,a path for pymol to readin input structure.
        '''
        self._set_structure(struct,filepath)


    def _set_reference_structure(self,struct:Union[str,Structure],filepath:Union[str,None]=None)->None:
        '''
        similar to _set_structure
        '''
        if isinstance(struct,str):
            self.reference=read_in(struct,struct)
            self.reference_path=struct
        elif isinstance(struct,Structure):
            assert isinstance(filepath,str),'must give me the path of pdb!'
            self.reference=struct.copy()
            self.reference_path=filepath
        else:
            raise TypeError
        dirty_seqs=sequence_from_object(self.reference,self.resname_dict)
        self.reference_sequences=wash_sequence(dirty_seqs)
    def _find_ligands(self):
        '''
        set the self.ligands, a list (may transfer to a Tuple later) containing all ligand Residue instance in ref.
        solvent will be removed by built-in rules of biopython.
        '''
        self.ligands:List[Residue]=[i for i in integrated_residue_iterator(self.reference) 
                        if i.id[0] not in [' ','W']]
    def _find_proteins(self):
        '''
        set the self.protein,a list (may transfer to a Tuple later) containing all protein Residue instance in ref
        '''
        self.proteins:List[Residue]=[i for i in integrated_residue_iterator(self.reference) 
                        if i.id[0]==' ']
    def _compare_chain(self):
        '''
        set self.chain_distance, a Dict[str,Dict[str,str]] object.
        outer key: chain id in object,inner key:chain id in reference
        value:distance calculate in dnd.
        '''
        self.chain_distance:Dict[str,Dict[str,str]]={}
        def process_dnd(dndvalue:str)->float:
            dnd_list=dndvalue.strip(';\n').strip('(').strip(')').split(',')
            # dnd_dict={i.split(':')[0]:i.split(':')[1] for i in dnd_list}
            distance=dnd_list[1].split(':')[1]
            # print(dndvalue,'\t',distance)
            return float(distance)

        with tempfile.TemporaryDirectory() as dirname:
            fasta_file=dirname+'/align.fasta'
            aln_file=dirname+'/align.aln'
            dnd_file=dirname+'/align.dnd'
            distance_dict={}
            for object_chain in self.sequences.keys():
                sub_dict={}
                for reference_chain in self.reference_sequences.keys():
                    write_fasta(prefix=self.object.id+'_',seq_dict=self.sequences,chainid=object_chain,output=fasta_file,rewrite=True)
                    write_fasta(prefix=self.reference.id+'_',seq_dict=self.reference_sequences,chainid=reference_chain,output=fasta_file)
                    clustalw_cline=ClustalwCommandline("clustalw2", infile=fasta_file,outfile=aln_file)
                    clustalw_cline()
                    with open(dnd_file,'r') as dnd:
                        dndvalue=''.join(dnd.readlines())
                    distance=process_dnd(dndvalue)
                    sub_dict[reference_chain]=distance
                distance_dict[object_chain]=sub_dict
        self.chain_distance=distance_dict

    def _align_sequence(self,object_chain:str,reference_chain:str)->pd.DataFrame:
        '''
        align two sequences, one from object structure and another from reference structure
        
        1.set self.alignment._records, a MultipleSeqAlignment instance.
        2.set self.alignment_frame, a DataFrame with 2 columns,  'object_residue' & 'reference_residue',
            each row is a residue mapping.
            when a residue is mapped to '-', util.RESIDUE_HOLDER object will be used to fill the blanket.

        Note:
        need to run self.map_neighbor to set self.mapped_neighbor under most cases!
        '''
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

    def align_sequence(self,object_chain:str,reference_chain:str)->pd.DataFrame:
        '''
        align two sequences, one from object structure and another from reference structure
        
        1.set self.alignment._records, a MultipleSeqAlignment instance.
        2.set self.alignment_frame, a DataFrame with 2 columns,  'object_residue' & 'reference_residue',
            each row is a residue mapping.
            when a residue is mapped to '-', util.RESIDUE_HOLDER object will be used to fill the blanket.

        Note:
        need to run self.map_neighbor to set self.mapped_neighbor under most cases!
        '''
        return self._align_sequence(object_chain,reference_chain)

    def _align_structure(self,save_ligand:bool=True,save_pocket:bool=False,
                            ligandfile:str='aligned_ligand.pdb',ref_pocketfile:str='aligned_pocket.pdb',obj_pocketfile:str='model_pocket.pdb')->Tuple[float,float]:
        '''
        invoke pymol to align the structure.
        load file of self.reference_path & self.object_path;
        select the ligand, ref_pocket & target_pocket, align the pockets
        (now use pymol.cmd.align but it's said that in low homology system, cmd.cealign/super performs better. )
        then save the ligand & pocket pdb upon request. 
 
        '''
        ref_model=self.reference_path.split('/')[-1].strip('.pdb')
        object_model=self.object_path.split('/')[-1].strip('.pdb')
        ligand_code=ref_model+' and '+reslist_to_guicode(self.selected_ligand,'pymol')
        ref_code=ref_model+' and '+reslist_to_guicode(self.neighbor,'pymol')
        model_code=object_model+' and '+reslist_to_guicode(self.mapped_neighbor,'pymol')
        pymol.cmd.load(self.object_path)
        pymol.cmd.load(self.reference_path)
        all_rmsd=pymol.cmd.align(mobile=ref_model,target=object_model)[0]
        pymol.cmd.select(name='lig',selection=ligand_code)
        pymol.cmd.select(name='ref_pocket',selection=ref_code)
        pymol.cmd.select(name='target_pocket',selection=model_code)
        _=pymol.cmd.align(mobile='ref_pocket',target='target_pocket')
        # print(_)
        pocket_rmsd=pymol.cmd.align(mobile='ref_pocket',target='target_pocket')[0]
        if save_ligand:
            pymol.cmd.save(ligandfile,selection='lig')
        if save_pocket:
            pymol.cmd.save(ref_pocketfile,selection='ref_pocket')
            pymol.cmd.save(obj_pocketfile,selection='target_pocket')
        pymol.cmd.delete('all')
        return all_rmsd,pocket_rmsd

    def align_structure(self,save_ligand:bool=True,save_pocket:bool=False,
                            ligandfile:str='aligned_ligand.pdb',ref_pocketfile:str='aligned_pocket.pdb',obj_pocketfile:str='model_pocket.pdb')->Tuple[float,float]:
        '''
        invoke pymol to align the structure.
        load file of self.reference_path & self.object_path;
        select the ligand, ref_pocket & target_pocket, align the pockets
        (now use pymol.cmd.align but it's said that in low homology system, cmd.cealign/super performs better. )
        then save the ligand & pocket pdb upon request. 
        '''
        return self._align_structure(save_ligand,save_pocket,
                            ligandfile,ref_pocketfile,obj_pocketfile)

    def init_mapped_neighbor(self):
        '''
        a tool to init and reset the `mapped_neighbor` 
        use when you mapped one system ,and want to run the next one.
        '''
        self.mapped_neighbor:List[Residue]=[]

    def map_neighbor(self):
        '''
        set self.mapped_neighbor, list of residues homologous with the self.neighbor in self.object.
        according to the alignment information in self.alignment_frame (gen by align_sequence),
        '''
        f=self.alignment_frame
        if not hasattr(self,'mapped_neighbor'):
            self.init_mapped_neighbor()
        else:
            self.mapped_neighbor.extend(f[f['reference_residue'].isin(self.neighbor)]['object_residue'].to_list())

    def load_reference(self,struct:Union[str,Structure]):
        '''
        set the self.reference & self.reference_sequences according to input structure.
        set self.ligands, a list contains all non-water ligand  residue in the reference.
        set self.proteins, a list contains all protein residue in the reference.
        set self.chain_distance, a Dict[str,Dict[str,str]] object.
        outer key: chain id in object,inner key:chain id in reference
        value:distance calculate in dnd.
        '''
        self._set_reference_structure(struct)
        self._find_ligands()
        self._find_proteins()
        # self._compare_chain()

    def set_ligand(self,ligand:Residue):
        '''
        choose a ligand from self.ligands
        set it to self.selected_ligand
        calculate protein residues within 6 A of self.selected_ligand
        set the self.neighbor as a list of these residues.
        set the self.neighbor_chain as a list of chain contains self.neighbor.
        '''
        assert ligand in self.ligands,'invalid ligand!'
        self.selected_ligand=ligand
        self.neighbor=residue_within_threshold(self.proteins,self.selected_ligand,6)
        self.neighbor_chain= tuple(set([i.get_parent().id for i in self.neighbor]))

    # def _produce_feature(self,object_chain:str,reference_chain:str):
    #     self._align_sequence(object_chain,reference_chain)

    def set_chain_as_ligand(self,chain:Chain):
        assert chain in list(self.reference.get_chains()),'invalid chain'
        self.selected_ligand=chain
        _other_residue=[ residue for residue in self.proteins if residue not in list(chain.get_residues()) ]
        self.neighbor=residue_within_threshold(_other_residue,self.selected_ligand,6)
        self.neighbor_chain= tuple(set([i.get_parent().id for i in self.neighbor]))

    def transform(self) -> Tuple[Tuple[str],Tuple[float]]:
        '''
        run transform when you are not sure about chain mapping relation ship.
        it will traverse all the combination, and return the pocket with the lowest pocket rmsd.
        chain map will be save in self.mapped_chain (tuple)
        best pocket rmsd will be save in self.rmsd
        the rmsd of whole structure and the pocket will be returned.
        '''
        object_chain=[i.id for i in self.object.get_chains()]
        num_pocket_chain=len(self.neighbor_chain)
        best_rmsd=(99999,99999)
        mapped_count=0
        best_combine=('place_holder',)
        for chain_combine in permutations(object_chain,num_pocket_chain):
            print(chain_combine,end='\t')
            self.init_mapped_neighbor()
            for i in range(num_pocket_chain):
                self._align_sequence(chain_combine[i],self.neighbor_chain[i])
                self.map_neighbor()
            
            _mapped_count=len([i for i in self.mapped_neighbor if i.resname!='XXX'])
            print(_mapped_count,end='\t')

            if _mapped_count>0:
                rmsd=self._align_structure(save_ligand=False,save_pocket=False)
            else:
                rmsd=(99999,99999)

            if rmsd[1]<best_rmsd[1] or (rmsd[1]<=best_rmsd[1]+0.01 and mapped_count<_mapped_count):
                best_rmsd=rmsd
                best_combine=chain_combine
                mapped_count=_mapped_count
            
            print(rmsd,end='\n')
            # print(best_rmsd)
        self.mapped_count=mapped_count
        self.mapped_chain=best_combine
        self.rmsd=best_rmsd
        return best_combine,best_rmsd

    def save_align(self,save_ligand:bool=True,save_pocket:bool=True,
                    ligandfile:str='aligned_ligand.pdb',ref_pocketfile:str='aligned_pocket.pdb',obj_pocketfile:str='model_pocket.pdb')->Tuple[float,float]:
        '''
        2 usages:
        1. if you have run self.tranform(), use this method to run alignment and save .pdb file
        2. if you are sure about the chain mapping relationship:
            1st, run `self.mapped_chain=tuple(str,str,...)` (object chain homological to chain in self.neighbor_chain)
            2nd, run this method to align (and save the structure, if you like).
        '''
        self.init_mapped_neighbor()
        for obj_chain,ref_chain in zip(self.mapped_chain,self.neighbor_chain):
            self._align_sequence(obj_chain,ref_chain)
            self.map_neighbor()
            self.mapped_count=len([i for i in self.mapped_neighbor if i.resname!='XXX'])
        return self._align_structure(save_ligand=save_ligand,save_pocket=save_pocket,
                                ligandfile=ligandfile,ref_pocketfile=ref_pocketfile,obj_pocketfile=obj_pocketfile)
        # raise NotImplementedError
        
        
        
