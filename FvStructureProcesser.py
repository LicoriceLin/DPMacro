import typing,warnings
# import numpy as np
import pandas as pd
import Bio.PDB as BP
import anarci
import re
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.Data import SCOPData
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio import BiopythonParserWarning,BiopythonWarning
allowd_scheme=typing.Literal['c','chothia','k','kabat','i','imgt','a','aho','m','martin','w','wolfguy']


class FvSpliter:
    '''
    '''

    def __init__(self)->None:
        """
        """
        pass

    def set_chain(self,pdb_object:Chain,include_hetatm:bool=True)->None:
        """
        """
        self.object=pdb_object.copy()
        self.object.detach_parent()
        # self.hetatm=Chain('hetatm')
        # if include_hetatm:
        #     for i in self.object:
        #         if i.id[0]!=' ':
        #             self.hetatm.add(i)
        #     for i in self.hetatm:
        #         self.object.detach_child(i.id)
        self.hetatm=FvSpliter._extract_hetatm(self.object)
        self.residue_map=pd.DataFrame()
        self.residue_map['original_id']=[residue.id for residue in self.object]
        self.residue_map.set_index(['original_id'],inplace=True)
        self.residue_map['residue_type']=[SCOPData.protein_letters_3to1[residue.get_resname()] 
                                            for residue in self.object]
        self.sequence=''.join(self.residue_map['residue_type'])
   
    def _run_anarci(self,scheme:allowd_scheme)->None:
        """
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

    @staticmethod
    def _extract_hetatm(input:Chain)->Chain:
        '''
        '''
        output=Chain('hetatm')
        for i in input:
            if i.id[0]!=' ':
                output.add(i)
        for i in output:
            input.detach_child(i.id)
        return output


class FvStructureProcesser:
    '''
    '''

    def __init__(self) -> None:
        '''
        '''
        pass
        
    def set_structure(self,file:str,id:str='tmp',scheme='a')->None:
        '''
        '''
        self.object=BP.PDBParser(QUIET=True).get_structure(id,file)
        self.scheme=scheme
        if len(self.object)>1:
            warnings.warn("multi-frame object detected."
             "only frame 0 will be processed.",
             BiopythonParserWarning,)
        self.chain_list={chain.id :''.join([SCOPData.protein_letters_3to1[residue.get_resname()] 
                                            for residue in chain 
                                            if residue.id[0]==' ']) 
                                            for chain in self.object[0]}
        
        self.Fv_count={chain:anarci.run_anarci(sequence,scheme=scheme)[1][0]
                        for chain,sequence in self.chain_list.items()}
        self.Fv_count={chain:len(count) if count is not None else 0
                        for chain,count in self.Fv_count.items()}

    def split_chains(self)->None:
        '''
        '''
        self.splited_chains={}
        for chain_id in self.chain_list:
            if self.Fv_count[chain_id]>0:
                spliter=FvSpliter()
                spliter.set_chain(self.object[0][chain_id])
                self.splited_chains[chain_id]=spliter.split(self.scheme)
            else:
                hetatm_chain=FvSpliter._extract_hetatm(self.object[0][chain_id])
                self.splited_chains[chain_id]={'Ag0':self.object[0][chain_id],'hetatm':hetatm_chain}
    
    def build_output(self,
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
                    FvStructureProcesser._add_chain_to_builder(
                            origin_chain['Ag0'],origin_chain_id,builder.structure[0])
                    if contain_hetatm:
                        FvStructureProcesser._add_chain_to_builder(
                                origin_chain['hetatm'],new_chain_id,builder.structure[0])
                elif Fv_count==1:
                    FvStructureProcesser._add_chain_to_builder(
                            origin_chain['Fv0'],origin_chain_id,builder.structure[0])
                    if contain_hetatm:
                        FvStructureProcesser._add_chain_to_builder(
                            origin_chain['Fv0'],new_chain_id,builder.structure[0])
                elif Fv_count==2:
                    FvStructureProcesser._add_chain_to_builder(
                            origin_chain['Fv0'],origin_chain_id,builder.structure[0])
                    FvStructureProcesser._add_chain_to_builder(
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
                            'check your system or split them into multiple file for porcess.',
                            BiopythonWarning,)
                            break
                        FvStructureProcesser._add_chain_to_builder(
                            segment,chain_ids[i],builder.structure[0])
                        i += 1
                if i==len(chain_ids):
                    break
        self.build_structure=builder.get_structure()
        return builder.get_structure()

    def write_processed_structure(self,file:str='tmp.pdb',write_end:bool=True, preserve_atom_numbering:bool=False)->None:
        '''
        '''
        io = BP.PDBIO()
        io.set_structure(self.build_structure)
        io.save(file,write_end=write_end,preserve_atom_numbering=preserve_atom_numbering)

    @staticmethod
    def _add_chain_to_builder(segment:Chain,new_id:str,structure:Structure)->None:
        '''
        '''
        segment_to_build=segment.copy()
        segment_to_build.id=new_id
        segment_to_build.detach_parent()
        structure.add(segment_to_build)

#test code
if __name__ == '__main__':
    import sys
    test_path=sys.argv[1]
    Fv=FvStructureProcesser()
    Fv.set_structure(test_path)
    Fv.split_chains()
    Fv.build_output(contain_linker=False,contain_antigen=False,contain_hetatm=False,reserve_chain_id=False)
    Fv.write_processed_structure('only_Fv_changed_id.pdb')
    Fv.build_output(reserve_chain_id=False)
    Fv.write_processed_structure('only_Fv.pdb')
    Fv.build_output(contain_linker=False,contain_antigen=False,contain_hetatm=False,reserve_chain_id=True)
    Fv.write_processed_structure('changed_id.pdb')
    Fv.build_output(reserve_chain_id=True)
    Fv.write_processed_structure('full.pdb')