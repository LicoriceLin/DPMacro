'''
some useful recurrent codes

'''


# import imp
import Bio.PDB as BP
import pandas as pd 
import re,os
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Entity import Entity
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import numpy as np
from .Data import amino3to1dict
# from distance_util import idstr2tuple
from typing import List,Tuple,Union,Generator,Iterable,Literal,Dict,Any,Callable
from collections.abc import Iterable as collections_Iterable

# allowd_scheme=Literal['c','chothia','k','kabat','i','imgt','a','aho','m','martin','w','wolfguy']
allowd_scheme=Literal['c','chothia','k','kabat','i','imgt','a','aho']
allowed_residue_source=Union[Entity,Iterable[Entity],Atom,Iterable[Atom]]
STRUCTURE_HOLDER=Structure('place_holder')
MODEL_HOLDER=Model('place_holder')
CHAIN_HOLDER=Chain('place_holder')
RESIDUE_HOLDER=Residue(('place_holder',0,' '),'XXX','')
ATOM_HOLDER=Atom('CA',np.array([0,0,0],dtype=np.float32),0.0,1.0,' ',' CA ',0,'C',None,None)
#set_parent has some unknown bug. use `Parent.add` instead. 
STRUCTURE_HOLDER.add(MODEL_HOLDER)
MODEL_HOLDER.add(CHAIN_HOLDER)
CHAIN_HOLDER.add(RESIDUE_HOLDER)
RESIDUE_HOLDER.add(ATOM_HOLDER)

def read_in(file:str,id:Union[str,None]=None)->Structure:
    '''
    read in a pdb file,return a Structure object
    id:the id of output Structure. Use path str as the Structure.id if id=None

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
    will be compatible with all `allowed_residue_source` in the future
    '''
    io = BP.PDBIO()
    io.set_structure(strcture)
    io.save(file,write_end=write_end,preserve_atom_numbering=preserve_atom_numbering)

def split_frame(file:str)->None:
    s=read_in(file)
    # if len(s)<=1:
    # 单个结构的？
    if len(s)<1:
        raise ValueError
    outdir=file.replace('.pdb','')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for i,frame in enumerate(s):
        write_out(frame,os.path.join(outdir,f'{i}.pdb'))

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

def impute_beta(object:allowed_residue_source,func:Callable[[Residue], float]):
    '''
    func needs a Residue as input, a float as output
    '''
    for atom in integrated_atom_iterator(object):
        atom.bfactor=func(atom.get_parent())

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
    a standard method to convert str/int to "standard triplet biopython residue code" tuple.
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
        
def impute_default_value(object:allowed_residue_source,key:str,default_value:Any)->None:
    '''
    used when not all Residue in `object` has a key of `key` in its .xtra dict.
    impute the given default value into this entry. 
    '''
    for residue in integrated_residue_iterator(object):
        residue.xtra[key]=residue.xtra.get(key,default_value)


def sequence_from_object(object:Union[Structure,Model],map:dict=amino3to1dict)->Dict[str,str]:
    '''
    only suggested run on Structure & Model instance
    return a {chainid:aa_seq} dict by map residue's `resname` property in map dict

    '''
    seq_dict=dict()
    for chain in object.get_chains():
        seq_dict[chain.id]=''.join([map.get(i.get_resname(),'X') for i in chain.get_residues()])
    return seq_dict

def sequence_from_frame(frame:pd.DataFrame,map:dict=amino3to1dict)->Dict[str,str]:
    '''
    must run on ResidueFeatureExtractor.frame properety

    '''
    def list3_to_str1(list3:Iterable[str])->str:
            return ''.join([map.get(i,'X') for i in list3])
    _seq_series=frame.groupby('chain')['resname'].apply(list3_to_str1)
    return dict(_seq_series)

def wash_sequence(seq_dict:Dict[str,str])->Dict[str,str]:
    '''
    use output_dict of `sequence_from_object` or `sequence_from_frame` as input.
    remove the 'X's at the beginning and end of all sequence.
    '''
    pure_dict={}
    for chain,seq in seq_dict.items():
        pure_dict[chain]=re.sub('X+$|^X+','',seq)
    return pure_dict

def integrated_residue_iterator(object:allowed_residue_source)->Generator[Residue,None,None]:
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
            yield from integrated_residue_iterator(i)
    else:
        raise TypeError

def split_multiframe(infile:str)->None:
    s=read_in(infile)
    i=0
    for model in s.get_models():
        write_out(model,infile.replace('.pdb',f'_{i}.pdb'))
        i +=1


def write_fasta(seq_dict:dict,prefix:str,chainid:Union[str,None]=None,output:str='poc.fasta',rewrite:bool=False):
    '''
    write out a fasta file from the seq_dict 
            ({chainid:seq},usually generated from util.sequence_from_object or util.sequence_from_frame)
    prefix:a str to write before chainid in seq_dict
    chainid which chain to write out. set None if you want to write all. default = None.
    output: file name of output. default = 'poc.fasta'
    rewrite: if clear the origin component of the output file. default = False.
    '''
    mode='w' if rewrite else 'a' 
    with open(output,mode) as f:
        if chainid is None:
            for chain,seq in seq_dict.items():
                f.write('>'+prefix+chain+'\n')
                f.write(seq+'\n')
        else:
            f.write('>'+prefix+chainid+'\n')
            f.write(seq_dict[chainid]+'\n')

def integrated_atom_iterator(object:Union[Entity,Iterable[Entity],Atom,Iterable[Atom]])->Generator[Atom,None,None]:
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
            yield from integrated_atom_iterator(i)
    else:
        raise TypeError

def _list_feature_into_residue(feature_list:Iterable,key:str,object:allowed_residue_source):
    '''
    suggest to use in some feature imputer.
    when the previous operation would produce a feature list(must be list),
    use this function to put these feature in object's Residues.
    '''
    i=0
    for residue in integrated_residue_iterator(object):
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

def renum_according_ref(infile:str,reffile:str,outfile:str)->None:
    '''
    '''
    s=read_in(infile)
    s_ref=read_in(reffile)
    assert len(list(integrated_residue_iterator(s)))==len(list(integrated_residue_iterator(s_ref))),'residue numbers are not consistent!'
    def tmp_id(x:tuple):
        return x[0],x[1],x[2]+'___'

    def true_id(x:tuple):
        return x[0],x[1],x[2].replace('___','')

    for res_,res_ref in zip(integrated_residue_iterator(s),integrated_residue_iterator(s_ref)):
        res_.id=tmp_id(res_ref.id)

    for res_ in integrated_residue_iterator(s):
        res_.id=true_id(res_.id)

    for chain_,chain_ref in zip(s.get_chains(),s_ref.get_chains()):
        chain_.id=chain_ref.id+'___'

    for chain_ in s.get_chains():
        chain_.id=chain_.id.replace('___','')

    write_out(s,outfile)

#test code
#need to be update
if __name__ == '__main__':
    import sys
    test_path=sys.argv[1]
    chain_for_extraction=sys.argv[2]
    structure=read_in(test_path)
    output=extract_hetatm(structure[0][chain_for_extraction],inplace=True)
    add_chain(output[1],'Z',structure[0])
    write_out(structure)