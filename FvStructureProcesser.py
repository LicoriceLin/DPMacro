'''
Objects to split and renumber residue/chains in Fv containing pdb.
available function: see the test codes.
'''
import typing,warnings
# import numpy as np
import pandas as pd
import Bio.PDB as BP
import anarci
import re
from util import read_in,write_out,extract_hetatm,add_chain

from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.Data import SCOPData
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio import BiopythonParserWarning,BiopythonWarning
allowd_scheme=typing.Literal['c','chothia','k','kabat','i','imgt','a','aho','m','martin','w','wolfguy']


class FvSpliter:
    '''
    semi-private object.
    a class only for Fv-containing chain process.
    useful method:
    .set_chain(pdb_object,include_hetatm) & .split(scheme)
    function obvious as the name indicated.  
    '''

    def __init__(self)->None:
        """
        """
        pass

    def set_chain(self,pdb_object:Chain,include_hetatm:bool=True)->None:
        """
        include_hetatm: if you expect chain contain hetatm residues(default=True). 
                        set to False when you want to omit these residue.
                        or when you want 
        """
        self.object=pdb_object.copy()
        self.object.detach_parent()
        self.hetatm=extract_hetatm(self.object,inplace=True)[1] if include_hetatm else Chain('hetatm')
        self.residue_map=pd.DataFrame()
        self.residue_map['original_id']=[residue.id for residue in self.object]
        self.residue_map.set_index(['original_id'],inplace=True)
        self.residue_map['residue_type']=[SCOPData.protein_letters_3to1[residue.get_resname()] 
                                            for residue in self.object]
        self.sequence=''.join(self.residue_map['residue_type'])
   
    def _run_anarci(self,scheme:allowd_scheme)->None:
        """
        private method.

        """
        self.scheme=scheme
        anarci_maps=anarci.run_anarci(self.sequence,scheme=scheme)[1][0]
        assert anarci_maps,'this chain contains no Fv fragment'
        id,end,anarci_order,Fv_annotation=0,-1,[],[]
        for anarci_map in anarci_maps:
            last_end,start,end=end,anarci_map[1],anarci_map[2]

            anarci_order.extend(
                [(' ',-1,' ')]*(start-last_end-1))
            Fv_annotation.extend(
                [f'Lk{id}']*(start-last_end-1))

            anarci_order.extend(
                [(' ',)+i[0] for i in anarci_map[0] if i[1]!='-'])
            Fv_annotation.extend(
                [f'Fv{id}']*(end-start+1))

            id += 1

        last_end,start=end,len(self.sequence)
        anarci_order.extend(
                [(' ',-1,' ')]*(start-last_end-1))
        Fv_annotation.extend(
                [f'Lk{id}']*(start-last_end-1))    
        self.residue_map['new_id']=anarci_order
        self.residue_map['annotation']=Fv_annotation

    def split(self,scheme:allowd_scheme)->typing.Dict[str,Chain]:
        '''
        '''
        self._run_anarci(scheme)
        self.splited_object={i:Chain(i) for i in self.residue_map['annotation'].unique()}
        self.splited_object['hetatm']=Chain('hetatm')
        for i in self.residue_map.iterrows():
            original_id,new_id,annotation=i[0],i[1]['new_id'],i[1]['annotation']
            self.object[original_id].segid=annotation
            residue=self.object[original_id].copy()
            residue.detach_parent()
            if 'Fv' in annotation:
                residue.id=new_id
            self.splited_object[annotation].add(residue)
        self.splited_object['hetatm']=self.hetatm
            #     self.split_object[i[1]['annotation']].append(self.object[i[0]].copy())
        return self.splited_object


class FvStructureProcesser:
    '''
    a class for whole structure containing Fv domain.
    not supporting multiframe right now. in this case you can use `FvSpliter`.
    '''

    def __init__(self,struct:typing.Union[str,Structure]) -> None:
        '''
        init the object,
        put the struct into .object property
        '''
        if isinstance(struct,str):
            self.object=read_in(struct,struct)
        elif isinstance(struct,Structure):
            self.object=struct.copy()
        else:
            raise TypeError
        
    def parse_structure(self,scheme='a')->None:
        '''
        scheme: the numbering scheme to split the chains.

        yielding following property:
        chain_list:dict,{chainid,sequence}
        Fv_count:dict,{chainid:how many Fv domain in each chain}
        splited_chains:{chainid:list[splited Chains object]} (from FvSpliter)
        '''
        self.scheme=scheme
        if len(self.object)>1:
            warnings.warn("multi-frame object detected."
             "only frame 0 will be processed.",
             BiopythonParserWarning,)
        self.chain_list={chain.id :''.join([SCOPData.protein_letters_3to1[residue.get_resname()] 
                                            for residue in chain 
                                            if residue.id[0]==' ']) 
                                            for chain in self.object[0] if len(chain)>0}
        
        self.Fv_count={chain:anarci.run_anarci(sequence,scheme=scheme)[1][0] if len(sequence)>0 else None
                        for chain,sequence in self.chain_list.items()}
        self.Fv_count={chain:len(count) if count is not None else 0
                        for chain,count in self.Fv_count.items()}
        self._split_chains()

    def _split_chains(self)->None:
        '''
        private method
        generate the 
        '''
        self.splited_chains={}
        for chain_id in self.chain_list:
            if self.Fv_count[chain_id]>0:
                spliter=FvSpliter()
                spliter.set_chain(self.object[0][chain_id])
                self.splited_chains[chain_id]=spliter.split(self.scheme)
            else:
                hetatm_chain=extract_hetatm(self.object[0][chain_id])
                self.splited_chains[chain_id]={'Ag0':self.object[0][chain_id],'hetatm':hetatm_chain}
    
    def build_whole(self,
                    structure_id:str='tmp',
                    contain_linker:bool=True,
                    contain_antigen:bool=True,
                    contain_hetatm:bool=True,
                    reserve_chain_id:bool=False)->Structure:
        '''
        '''
        if reserve_chain_id and max(self.Fv_count.values())>2:
            warnings.warn('some chain contains more than 2 Fv fragments, '
                        'which means `reserve_chain_id` wont work. ',
                        BiopythonWarning,)
            reserve_chain_id=False
        if reserve_chain_id and contain_linker:
            warnings.warn('reserve chainid while contain linker will be too complicated but trivial, '
            'this function has not been implemented',
                        BiopythonWarning,)
            contain_linker=False

        id_match_string='(Fv'
        if contain_linker:
            id_match_string += '|Lk'
        if contain_antigen:
            id_match_string += '|Ag'
        if contain_hetatm:
            id_match_string += '|hetatm'
        id_match_string += ')'
        id_pattern=re.compile(id_match_string)

        builder=StructureBuilder()
        builder.init_structure(structure_id)
        builder.init_model(0)
        if reserve_chain_id:
            for origin_chain_id,origin_chain in self.splited_chains.items():
                new_chain_id=origin_chain_id.upper() if origin_chain_id.islower() else origin_chain_id.lower()
                Fv_count=self.Fv_count[origin_chain_id]
                if Fv_count==0 and contain_antigen:
                    add_chain(
                            origin_chain['Ag0'],origin_chain_id,builder.structure[0])
                    if contain_hetatm:
                        add_chain(
                                origin_chain['hetatm'],new_chain_id,builder.structure[0])
                elif Fv_count==1:
                    add_chain(
                            origin_chain['Fv0'],origin_chain_id,builder.structure[0])
                    if contain_hetatm:
                        add_chain(
                            origin_chain['Fv0'],new_chain_id,builder.structure[0])
                elif Fv_count==2:
                    add_chain(
                            origin_chain['Fv0'],origin_chain_id,builder.structure[0])
                    add_chain(
                            origin_chain['Fv1'],new_chain_id,builder.structure[0])
                    if contain_hetatm:
                        warnings.warn('two Fv fragment has occupied all the chain id, ' 
                            'so the hetatm will be left out.',
                            BiopythonWarning,)
        else:
            chain_ids=[chr(i) for i in list(range(65,91))+list(range(97,123))+list(range(48,58)) ]     
            i=0   
            for origin_chain in self.splited_chains.values():
                for segment in origin_chain.values():
                    if id_pattern.match(segment.id):
                        if i==len(chain_ids):
                            warnings.warn('builder has run out of all available chain ids'
                            'check your system or turn to the more personalized build_single method.',
                            BiopythonWarning,)
                            break
                        add_chain(
                            segment,chain_ids[i],builder.structure[0])
                        i += 1
                if i==len(chain_ids):
                    break
        self.built_structure=builder.get_structure()
        return builder.get_structure()

    def build_single(self,
                        chainid:str,
                        allocate_chain_ids:list=[chr(i) for i in list(range(65,91))+list(range(97,123))],
                        structure_id:typing.Union[str,None]=None,
                        )->Structure:
        '''
        split
        chainid:chain to build
        allocate_chain_ids: a list contains the chain id you want to allocate to the splited chains.
        structure_id : the structure id for the built new structure. None to keep the original id.
        '''
        allocate_chain_ids=[i for i in allocate_chain_ids if (i not in self.Fv_count.keys()) or i==chainid ]
        assert len(self.splited_chains[chainid])<=len(allocate_chain_ids),'more chain id needed!'
        
        builder=StructureBuilder()
        structure_id= structure_id if structure_id else self.object.id
        builder.init_structure(structure_id)
        builder.init_model(0)
        
        for chain in self.object[0]:
            if chain.id != chainid:
                add_chain(chain,chain.id,builder.structure[0])
        
        i=0
        for segment in self.splited_chains[chainid].values():
            add_chain(segment,allocate_chain_ids[i],builder.structure[0])
            i += 1
        self.built_structure=builder.get_structure()
        return builder.get_structure()

    def update_structure(self)->None:
        '''
        '''
        self.object=self.built_structure
        self.parse_structure(self.scheme)
  
    def write_processed_structure(self,file:str='tmp.pdb',write_end:bool=True, preserve_atom_numbering:bool=False)->None:
        '''
        '''
        write_out(self.built_structure,file,write_end,preserve_atom_numbering)


#test code
if __name__ == '__main__':
    import sys,os
    #arg1:the path to test pdb file
    test_path=sys.argv[1]
    #arg2&3:the first&second chain to test in `build_single` method
    test_chain1=sys.argv[2]
    test_chain2=sys.argv[3]

    os.mkdir('test')
    Fv=FvStructureProcesser(test_path)
    Fv.parse_structure()
    
    #test function:write only Fvs, renumber chain id from A 
    Fv.build_whole(contain_linker=False,contain_antigen=False,contain_hetatm=False,reserve_chain_id=False)
    Fv.write_processed_structure('test/only_Fv_changed_id.pdb')
    
    #test function:write all component, renumber chain id from A 
    Fv.build_whole(reserve_chain_id=False)
    Fv.write_processed_structure('test/only_Fv.pdb')

    #test function:write only Fvs,keep the original chain id ( upper & lower style for diabody )
    Fv.build_whole(contain_linker=False,contain_antigen=False,contain_hetatm=False,reserve_chain_id=True)
    Fv.write_processed_structure('test/keep_id.pdb')

    #test function:write all component, keep the original chain id ( upper & lower style for diabody )
    #can be unstable in many cases.
    Fv.build_whole(reserve_chain_id=True)
    Fv.write_processed_structure('test/full.pdb')
    
    #test function:split chain1, renumber from chain A ,and skip the already existing chains
    Fv.build_single(test_chain1)
    Fv.write_processed_structure('test/single_chain.pdb')
    
    #test function:if update the Fv.object according Fv.built_structure runs well. 
    Fv.update_structure()
    Fv.build_single(test_chain2)
    Fv.write_processed_structure('test/iterated.pdb')