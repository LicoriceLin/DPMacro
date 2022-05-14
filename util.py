'''
some useful stand-alone tools

cmds for condes validation:
python utils arg1 arg2
arg1 : path of the pdb testcase
arg2 : the chain to run extract_hetatm on 

the output will be like: 
the arg2 chain's hetatms will be renamed to chain Z and put behind all other chains.

more function will to be added will the progression of reconstruction.

'''


import Bio.PDB as BP
import pandas as pd 
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Entity import Entity
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from Data import amino3to1dict
# from distance_util import idstr2tuple
from typing import List,Tuple,Union,Generator,Iterable,Literal,Dict,Any
from collections.abc import Iterable as collections_Iterable

# allowd_scheme=Literal['c','chothia','k','kabat','i','imgt','a','aho','m','martin','w','wolfguy']
allowd_scheme=Literal['c','chothia','k','kabat','i','imgt','a','aho']
allowed_residue_source=Union[Entity,Iterable[Entity],Atom,Iterable[Atom]]
STRUCTURE_HOLDER=Structure('place_holder')
MODEL_HOLDER=Model('place_holder')
CHAIN_HOLDER=Chain('place_holder')
RESIDUE_HOLDER=Residue(('place_holder',0,' '),'XXX','')
MODEL_HOLDER.set_parent(STRUCTURE_HOLDER)
CHAIN_HOLDER.set_parent(MODEL_HOLDER)
RESIDUE_HOLDER.set_parent(CHAIN_HOLDER)

def read_in(file:str,id:Union[str,None]=None)->Structure:
    '''
    Parse a PDB file into a Structure object in `PDBParser`'s default mode
    allow an extra parameter `id` for the structure's name. 
    '''
    if isinstance(id,str):
        return BP.PDBParser(QUIET=True).get_structure(id,file)
    else:
        return BP.PDBParser(QUIET=True).get_structure(file,file)

def write_out(strcture:Entity,file:str='tmp.pdb',write_end:bool=True, preserve_atom_numbering:bool=False)->None:
    '''
    write out an Entity (from the whole structure to a single atom) to pdb file.
    use the `write_end` to write an END line, 
    the `preserve_atom_numbering` to renumber atom from 1.
    '''
    io = BP.PDBIO()
    io.set_structure(strcture)
    io.save(file,write_end=write_end,preserve_atom_numbering=preserve_atom_numbering)

def extract_hetatm(input:Chain,inplace:bool=False)->List[Chain]:  
    '''
    split one chain into two, one contains only ATOM lines, the other contains only HETATM lines.
    allow an extra parameter `inplace`. hetatms in original chain will be removed if it is true.
    '''
    output:List[Chain]=[]
    output.append(Chain('atom'))
    output.append(Chain('hetatm'))
    for i in input:
        i_copy=i.copy()
        if i.id[0]==' ':
            output[0].add(i_copy)
        else:
            output[1].add(i_copy)
    if inplace:
        for i in output[1]:
            input.detach_child(i.id)
    return output

def add_chain(segment:Chain,new_id:str,Model:Model)->None:
    '''
    add an exist chain to a Model object (a `frame` in a structure)
    '''
    segment_to_build=segment.copy()
    segment_to_build.id=new_id
    segment_to_build.detach_parent()
    Model.add(segment_to_build)

def to_resid(input:Union[int,str,Tuple[str,int,str]])->Tuple[str,int,str]:
    '''
    a standard method to convert str, int to "standard triplet biopython residue code" tuple.
    return itself if the standard code is give.
    '''
    if isinstance(input,tuple) and len(input)==3:
        try:
            return tuple(str(input[0]),int(input[1]),str(input[2]))
        except:
            raise TypeError
    else:
        try:
            return (' ',int(input),' ')
        except:
            raise TypeError
        
def _impute_default_value(object:allowed_residue_source,key:str,default_value:Any)->None:
    for residue in _integrated_residue_iterator(object):
        residue.xtra[key]=residue.xtra.get(key,default_value)


def sequence_from_object(object:Union[Structure,Model],map:dict=amino3to1dict)->Dict[str,str]:
    '''
    only suggested run on Structure & Model instance
    return a {chainid:aa_seq} dict

    '''
    seq_dict=dict()
    for chain in object.get_chains():
        seq_dict[chain.id]=''.join([map.get(i.get_resname(),'X') for i in chain.get_residues()])
    return seq_dict

def sequence_from_frame(frame:pd.DataFrame,map:dict=amino3to1dict)->Dict[str,str]:
    '''
    must run on ResidueFeatureExtractor.object

    '''
    def list3_to_str1(list3:Iterable[str])->str:
            return ''.join([map.get(i,'X') for i in list3])
    _seq_series=frame.groupby('chain')['resname'].apply(list3_to_str1)
    return dict(_seq_series)

def _integrated_residue_iterator(object:allowed_residue_source)->Generator[Residue,None,None]:
    '''
    iterate the Residue instance in Entity, Atom(return its parent in this case ),
        AND anything iterable composed of these two kind of instance
    '''
    if isinstance(object,Entity) and object.level!='R':
        yield from object.get_residues()
    elif isinstance(object,Entity) and object.level=='R':
        yield from [object]
    elif isinstance(object,Atom):
        yield from [object.get_parent()]
    elif isinstance(object,collections_Iterable):
        for i in object:
            yield from _integrated_residue_iterator(i)
    else:
        raise TypeError

def _integrated_atom_iterator(object:Union[Entity,Iterable[Entity],Atom,Iterable[Atom]])->Generator[Atom,None,None]:
    '''
        iterate the Atom instance in Entity, Atom(return their parent in this case ),
        AND anything iterable composed of these two kind of instance
    '''
    if isinstance(object,Entity):
        yield from object.get_atoms()
    elif isinstance(object,Atom):
        yield from [object]
    elif isinstance(object,collections_Iterable):
        for i in object:
            yield from _integrated_atom_iterator(i)
    else:
        raise TypeError

def _list_feature_into_residue(feature_list:Iterable,key:str,object:allowed_residue_source):
    '''
    suggest to use in some feature imputer.
    when the previous operation would produce a feature list(must be list),
    use this function to put these feature in object's Residues.
    '''
    i=0
    for residue in _integrated_residue_iterator(object):
        residue.xtra[key]=feature_list[i]
        i+=1
    if i!=len(feature_list):
        raise ValueError
    
def _list_feature_into_frame(feature_list:Iterable,key:str,frame:pd.DataFrame):
    '''
    suggest to use in some feature imputer.
    when the previous operation would produce a feature list(must be list),
    use this function to put these feature into the frame 
    commonly the ResidueFeatureExtractor.object)
    '''
    if len(frame)!=len(feature_list):
        raise ValueError
    else:
        frame[key]=feature_list


#test code
if __name__ == '__main__':
    import sys
    test_path=sys.argv[1]
    chain_for_extraction=sys.argv[2]
    structure=read_in(test_path)
    output=extract_hetatm(structure[0][chain_for_extraction],inplace=True)
    add_chain(output[1],'Z',structure[0])
    write_out(structure)