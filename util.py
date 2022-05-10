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
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Entity import Entity
from Bio.PDB.Atom import Atom
import pandas as pd 
from Data import amino3to1dict
from typing import List,Tuple,Union

def read_in(file:str,id:str='tmp')->Structure:
    '''
    Parse a PDB file into a Structure object in `PDBParser`'s default mode
    allow an extra parameter `id` for the structure's name. 
    '''
    return BP.PDBParser(QUIET=True).get_structure(id,file)

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
    a standard method to convert str, int to standard triplet biopython residue code tuple.
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
        
def idtuple2str(idtuple:Union[Tuple,str,int])->str:
    '''
    '''
    if isinstance(idtuple,tuple):
        idlist=[str(i) if i != ' ' else '_'  for i in idtuple  ]
        return '|'.join(idlist)
    else:
        return str(idtuple) if idtuple != ' ' else '_'
    
def idstr2tuple(idstr:str)->Tuple:
    '''
    '''
    def _(str):
        try:
            out=int(str)
        except ValueError:
            if str=='_':
                out=' '
            else:
                out=str
        return out
    idlist=idstr.split('|')
    idlist=[_(i) for i in idlist]
    if len(idlist)>1:
        return tuple(idlist)
    elif len(idlist)==1:
        return idlist[0]
    else:
        return ValueError

def atom_id_from_frame(frame:pd.DataFrame,index:int)->Tuple[int,str,tuple,tuple]:
    '''
    '''
    assert set(['model','chain','residue','atom'])<set(frame.columns),'frame does not contain adequate columns'
    _id=frame.loc[index][['model','chain','residue','atom']]
    id=tuple(_id.map(idstr2tuple))
    return id

def atom_object_from_frame(frame:pd.DataFrame,index:int)->Atom:
    '''
    '''
    assert set(['object'])<set(frame.columns),'frame does not contain adequate columns'
    object=frame.loc[index]['object']
    return object

def residue_id_from_frame(frame:pd.DataFrame,index:int)->Tuple[int,str,tuple,tuple]:
    '''
    '''
    assert set(['model','chain','residue'])<set(frame.columns),'frame does not contain adequate columns'
    _id=frame.loc[index][['model','chain','residue']]
    id=tuple(_id.map(idstr2tuple))
    return id

def get_seq(struct:Structure):
    seq_dict=dict()
    for chain in struct.get_chains():
        seq_dict[chain.get_full_id()[1:]]=''.join([i.getresname() for i in chain.get_residues()])
    return seq_dict


#test code
if __name__ == '__main__':
    import sys
    test_path=sys.argv[1]
    chain_for_extraction=sys.argv[2]
    structure=read_in(test_path)
    output=extract_hetatm(structure[0][chain_for_extraction],inplace=True)
    add_chain(output[1],'Z',structure[0])
    write_out(structure)