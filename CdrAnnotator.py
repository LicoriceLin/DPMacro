'''
antibody annotation functions.
run anarci and get the new residue id codes,
annotate the CDRs and their vicinity.
'''

from ctypes import Structure
from typing import Dict, Tuple,Union,List,Iterable
# import warnings

import anarci
import pandas as pd
from itertools import chain
from numpy.linalg import norm

from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Polypeptide import PPBuilder

from .BaseClasses import ResidueFeatureExtractor
from .util import _list_feature_into_residue,_list_feature_into_frame,integrated_residue_iterator
from .util import write_out,add_chain
from .util import allowed_residue_source,allowd_scheme
from .Data import CDR_annotations
# import distance_util as du



def run_anarci(sequence_dict:Dict[str,str],scheme:allowd_scheme='a',
                object:Union[allowed_residue_source,None]=None,frame:Union[pd.DataFrame,None]=None)->Tuple[list,list,list]:
    """
    use sequence_dict ( generated by ResidueFeatureExtractor._produce_sequence ) rather than structures as input.
    

    return 3 list:
        anarci_order((' ',112,'A'),...);
        Fv_annotation(Fv1,Lk0,...);
        chain_type(H,L,K,X)

    you can give the object or frame as parameter if you want to impute the output lists into them.
    """
    anarci_order,Fv_annotation,chain_type=[],[],[]
    for seq in sequence_dict.values():
        assert 'X' not in seq,f'invalid seqs: {seq}'
        anarci_results=anarci.run_anarci(seq,scheme=scheme)
        anarci_maps=anarci_results[1][0]
        chain_infos=anarci_results[2][0]
        # assert anarci_maps,'this chain contains no Fv fragment'
        if anarci_maps:
            id,end,_anarci_order,_Fv_annotation,_chain_type=0,-1,[],[],[]
            for anarci_map,chain_info in zip(anarci_maps,chain_infos):
                last_end,start,end=end,anarci_map[1],anarci_map[2]

                _anarci_order.extend(
                    [(' ',-1,' ')]*(start-last_end-1))
                _Fv_annotation.extend(
                    [f'Lk{id}']*(start-last_end-1))
                _chain_type.extend(
                    ['X']*(start-last_end-1))

                _anarci_order.extend(
                    [(' ',)+i[0] for i in anarci_map[0] if i[1]!='-'])
                _Fv_annotation.extend(
                    [f'Fv{id}']*(end-start+1))
                _chain_type.extend(
                    [chain_info['chain_type']]*(end-start+1))    

                id += 1

            last_end,start=end,len(seq)
            _anarci_order.extend(
                    [(' ',-1,' ')]*(start-last_end-1))
            _Fv_annotation.extend(
                    [f'Lk{id}']*(start-last_end-1))
            _chain_type.extend(
                    ['X']*(start-last_end-1))    
            anarci_order.extend(_anarci_order)
            Fv_annotation.extend(_Fv_annotation)
            chain_type.extend(_chain_type)

        else:
            anarci_order.extend([(' ',-1,' ')]*len(seq))
            Fv_annotation.extend(['Lk0']*len(seq))
            chain_type.extend(['X']*len(seq))

        
    if isinstance(object,Entity):
        _list_feature_into_residue(anarci_order,'anarci_order',object)
        _list_feature_into_residue(Fv_annotation,'Fv_annotation',object)
        _list_feature_into_residue(chain_type,'chain_type',object)
    
    if isinstance(frame,pd.DataFrame):
        _list_feature_into_frame(anarci_order,'anarci_order',frame)
        _list_feature_into_frame(Fv_annotation,'Fv_annotation',frame)
        _list_feature_into_frame(chain_type,'chain_type',frame)

    return anarci_order,Fv_annotation,chain_type

def impute_cdr(object:allowed_residue_source,scheme:allowd_scheme='a')->None:
    '''
    must run `_run_anarci` before this function
    for structure which has been imputed 'anarci_order',
    impute `CDR` annotation (CDR3,FW1,AC(anchor)3,LK(linker),...) 
    '''
    def key_func(scheme:str,chaintype:str)->str:
        if scheme in ['a','i','aho','imgt']:
            return scheme
        elif scheme in ['c','chothia','k','kabat']:
            return scheme+'|'+chaintype if chaintype in ['H','L','K'] else 'a'
        else:
            raise ValueError
    for i in integrated_residue_iterator(object):
        annotation=CDR_annotations[key_func(scheme,i.xtra['chain_type'])]
        if i.xtra['anarci_order'][1]<annotation['FR1'] or i.xtra['anarci_order'][1]>annotation['end']:
            i.xtra['CDR']='LK'

        elif i.xtra['anarci_order'][1]>=annotation['FR1'] and i.xtra['anarci_order'][1]<annotation['CDR1']-1:
            i.xtra['CDR']='FW1'
        elif i.xtra['anarci_order'][1] in [annotation['CDR1']-1,annotation['FR2']] :
            i.xtra['CDR']='AC1'
        elif i.xtra['anarci_order'][1]>=annotation['CDR1'] and i.xtra['anarci_order'][1]<annotation['FR2']:
            i.xtra['CDR']='CDR1'

        elif i.xtra['anarci_order'][1]>annotation['FR2'] and i.xtra['anarci_order'][1]<annotation['CDR2']-1:
            i.xtra['CDR']='FW2'
        elif i.xtra['anarci_order'][1] in [annotation['CDR2']-1,annotation['FR3']] :
            i.xtra['CDR']='AC2'
        elif i.xtra['anarci_order'][1]>=annotation['CDR2'] and i.xtra['anarci_order'][1]<annotation['FR3']:
            i.xtra['CDR']='CDR2'

        elif i.xtra['anarci_order'][1]>annotation['FR3'] and i.xtra['anarci_order'][1]<annotation['CDR3']-1:
            i.xtra['CDR']='FW3'
        elif i.xtra['anarci_order'][1] in [annotation['CDR3']-1,annotation['FR4']] :
            i.xtra['CDR']='AC3'
        elif i.xtra['anarci_order'][1]>=annotation['CDR3'] and i.xtra['anarci_order'][1]<annotation['FR4']:
            i.xtra['CDR']='CDR3'

        elif i.xtra['anarci_order'][1]>annotation['FR4'] and i.xtra['anarci_order'][1]<=annotation['end']:
            i.xtra['CDR']='FW4'

def peptide_is_fv(list)->bool:
    '''
    temparory, not rubust function to distinguish fv from lk
    use peptide list as input
    '''
    if 90<len(list)<130:
        return True
    else:
        return False

def actually_same_chain(pep1:List[Residue],pep2:List[Residue])->bool:
    pep1_n=pep1[-1]['N']
    pep2_ca=pep2[0]['CA']
    dis=norm((pep1_n.coord-pep2_ca.coord),ord=2)
    if dis >5:
        return False
    else:
        return True

def consecutive_resid_list(start_id:int,len)->List[Tuple[str,int,str]]:
    resid_list=[]
    for i in range(start_id,start_id+len):
        resid_list.append((' ',i,' '))
    return resid_list

def renumber_from_resid_list(pep:List[Residue],resid_list:List[Tuple[str,int,str]])->List[Residue]:
    if len(pep) != len(resid_list):
        raise ValueError('check the lenth of `resid_list`')
    new_pep=[]
    for res,id in zip(pep,resid_list):
        new_res=res.copy()
        new_res.id=id
        new_pep.append(new_res)
    return new_pep

def build_new_chain_from_pep(peps:Iterable[List[Residue]],newchainid:str)->Chain:
    new_chain=Chain(newchainid)
    for residue in chain(*peps):
        # print(residue)
        res_=residue.copy()
        new_chain.add(res_)
    return new_chain

def lk_fv(pep1:List[Residue],pep2:List[Residue],id:str)->Chain:
    id_list=consecutive_resid_list(pep2[0].id[1]-len(pep1),len(pep1))
    chain_=build_new_chain_from_pep(
        [renumber_from_resid_list(pep1,id_list),pep2],
        id
    )
    return chain_

def fv_lk(pep1:List[Residue],pep2:List[Residue],id:str)->Chain:
    id_list=consecutive_resid_list(pep1[-1].id[1]+1,len(pep2))
    chain_=build_new_chain_from_pep(
        [pep1,renumber_from_resid_list(pep2,id_list)],
        id
    )
    return chain_

def lk_fv_lk(pep1:List[Residue],pep2:List[Residue],pep3:List[Residue],id:str)->Chain:
    id_list1=consecutive_resid_list(pep2[0].id[1]-len(pep1),len(pep1))
    id_list3=consecutive_resid_list(pep2[-1].id[1]+1,len(pep3))
    chain_=build_new_chain_from_pep(
        [   renumber_from_resid_list(pep1,id_list1),
            pep2,
            renumber_from_resid_list(pep3,id_list3)],
        id
    )
    return chain_

def fv(pep1:List[Residue],id:str)->Chain:
    return build_new_chain_from_pep([pep1],id)

def build_structure(chain_dict:Dict[str,Chain],id:str='tmp')->Structure:
    builder=StructureBuilder()
    builder.init_structure(structure_id=id)
    builder.init_model(0)

    for chainid,chain in chain_dict.items():
        add_chain(chain,chainid,builder.structure[0])
    
    return builder.get_structure()


class hmt_FvProcessor(ResidueFeatureExtractor):
    def __init__(self, scheme:allowd_scheme) -> None:
        self.scheme=scheme
        ResidueFeatureExtractor.__init__(self, 'FvProcessor', True)
    def _produce_feature(self) -> None:
        self._produce_sequence()
        run_anarci(self.sequences,self.scheme,self.object)
        impute_cdr(self.object,self.scheme)
        self._object_feature_to_frame()
    def renumber_structure(self)->Structure:
        last_Fv_annotation='holder'
        chain_dict={}
        
        def set_new_chain_id(new_chain_type:str,chain_dict:dict,linker_code:list=['A','B','C','D'])->str:
            if new_chain_type == 'H':
                new_chain_code='H'
            elif new_chain_type in ['K','L']:
                new_chain_code='L'
            elif new_chain_type=='X':
                i=0
                while linker_code[i] in chain_dict.keys():
                    i+=1
                new_chain_code=linker_code[i]
            return new_chain_code


        def set_bfactor(CDR:str,chaintype:str)->float:
            if chaintype in ['K','L']:
                if CDR=='CDR1':
                    return 14.25
                elif CDR=='CDR2':
                    return 28.5
                elif CDR=='CDR3':
                    return 42.75
                else:
                    return 0
            elif chaintype=='H':
                if CDR=='CDR1':
                    return 57
                elif CDR=='CDR2':
                    return 71.25
                elif CDR=='CDR3':
                    return 85.5
                else:
                    return 0
            else:
                return 0
        
        for residue in integrated_residue_iterator(self.object):
            #init new chain if the `Fv_annotation` is changed.
            if residue.xtra['Fv_annotation']!=last_Fv_annotation:
                last_Fv_annotation=residue.xtra['Fv_annotation']
                chain_id=set_new_chain_id(residue.xtra['chain_type'],chain_dict)
                chain_dict[chain_id]=Chain(chain_id)
                if chain_id not in ['H','L']:
                    linker_resid=1
            #set residue id
            residue_copy=residue.copy()
            if chain_id in ['H','L']:
                residue_copy.id=residue_copy.xtra['anarci_order']
            else:
                residue_copy.id=(' ',linker_resid,' ')
                linker_resid+=1

            #set bfactor
            # atom:Atom=Atom()
            for atom in residue_copy.get_atoms():
                atom.bfactor=set_bfactor(residue_copy.xtra['CDR'],residue_copy.xtra['chain_type'])
            
            chain_dict[chain_id].add(residue_copy)

        #build structure
        # builder=StructureBuilder()
        # builder.init_structure(structure_id='renumbered_Fv')
        # builder.init_model(0)

        # for chainid,chain in chain_dict.items():
        #     add_chain(chain,chainid,builder.structure[0])
        
        # self.built_structure=builder.get_structure()
        self.built_structure=build_structure(chain_dict,id='renumbered_Fv')
        return self.built_structure

    def recombine_chain(self)->Structure:
        if len(self.built_structure[0])==2:
            return self.built_structure[0]
        else:
            ppb=PPBuilder(radius=5)
            peps=ppb.build_peptides(self.built_structure)
            judge_fv=[peptide_is_fv(i) for i in peps]

            if judge_fv == [False,True,True]:
                chainH=lk_fv(peps[0],peps[1],'H')
                chainL=fv(peps[2],'L')

            elif judge_fv == [True,False,True]:
                if actually_same_chain(peps[0],peps[1]):
                    chainH=fv_lk(peps[0],peps[1],'H')
                    chainL=fv(peps[2],'L')
                elif actually_same_chain(peps[1],peps[2]):
                    chainH=fv(peps[0],'H')
                    chainL=lk_fv(peps[1],peps[2],'L')
                else:
                    raise ValueError

            elif judge_fv == [True,True,False]:
                chainH=fv(peps[0],'H')
                chainL=fv_lk(peps[1],peps[2],'L')

            elif judge_fv == [False, True, False, True]:
                if actually_same_chain(peps[1],peps[2]):
                    chainH=lk_fv_lk(peps[0],peps[1],peps[2],'H')
                    chainL=fv(peps[3],'L')
                elif actually_same_chain(peps[1],peps[2]):
                    chainH=lk_fv(peps[0],peps[1],'H')
                    chainL=lk_fv(peps[2],peps[3],'L')
                else:
                    raise ValueError

            elif judge_fv == [False, True, True, False]:
                chainH=lk_fv(peps[0],peps[1],'H')
                chainL=fv_lk(peps[2],peps[3],'L')

            elif judge_fv == [True, False, True, False]:
                if actually_same_chain(peps[0],peps[1]):
                    chainH=fv_lk(peps[0],peps[1],'H')
                    chainL=fv_lk(peps[2],peps[3],'L')
                elif actually_same_chain(peps[1],peps[2]):
                    chainH=fv(peps[0],'H')
                    chainL=lk_fv_lk(peps[1],peps[2],peps[3],'L')
                else:
                    raise ValueError

            elif judge_fv == [True, False, False, True]:
                chainH=fv_lk(peps[0],peps[1],'H')
                chainL=fv_lk(peps[2],peps[3],'L')
            elif judge_fv == [False, True, False, False, True]:
                chainH=lk_fv_lk(peps[0],peps[1],peps[2],'H')
                chainL=lk_fv(peps[3],peps[4],'L')
            elif judge_fv == [False, True, False, True, False ]:
                if actually_same_chain(peps[1],peps[2]):
                    chainH=lk_fv_lk(peps[0],peps[1],peps[2],'H')
                    chainL=fv_lk(peps[3],peps[4],'L')
                elif actually_same_chain(peps[2],peps[3]):
                    chainH=lk_fv(peps[0],peps[1],'H')
                    chainL=lk_fv_lk(peps[2],peps[3],peps[4],'L')
                else:
                    raise ValueError
            elif judge_fv == [True, False, False, True, False ]:
                chainH=fv_lk(peps[0],peps[1],'H')
                chainL=lk_fv_lk(peps[2],peps[3],peps[4],'L')
            elif judge_fv == [False, True, False, False, True, False]:
                chainH=lk_fv_lk(peps[0],peps[1],peps[2],'H')
                chainL=lk_fv_lk(peps[3],peps[4],peps[5],'L')
            else:
                raise ValueError
            struct=build_structure({'H':chainH,'L':chainL},'integrate_Fv')
            self.built_structure=struct
            return struct

    def write_built_structure(self,filename:str):
        write_out(self.built_structure,filename)

    def write_scheme_csv(self,filename:str):
        self.object=self.built_structure
        resid_list,resname_list,chainid_list,icode_list,cdr_list=[],[],[],[],[]
        
        for residue in integrated_residue_iterator(self.object):
            resid_list.append(residue.id[1])
            resname_list.append(residue.resname)
            chainid_list.append(residue.get_parent().id)
            icode_list.append(residue.id[2])
            cdr_list.append(residue.xtra['CDR'])

        output=pd.DataFrame()
        output['resid']=resid_list
        output['resname']=resname_list
        output['chainid']=chainid_list
        output['icode']=icode_list
        output['CDR']=cdr_list
        output.to_csv(filename,index=False)
        # self.frame=pd.DataFrame()
        # self._object_feature_to_frame()
        # self.frame['CDR'].to_json(filename)

            

