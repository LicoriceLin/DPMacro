import Bio.PDB as BP
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Entity import Entity
from typing import List

def read_in(file:str,id:str='tmp')->Structure:
    '''
    '''
    return BP.PDBParser(QUIET=True).get_structure(id,file)

def write_out(strcture:Entity,file:str='tmp.pdb',write_end:bool=True, preserve_atom_numbering:bool=False)->None:
    '''
    '''
    io = BP.PDBIO()
    io.set_structure(strcture)
    io.save(file,write_end=write_end,preserve_atom_numbering=preserve_atom_numbering)

def extract_hetatm(input:Chain,inplace=False)->List[Chain]:  
    '''
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

def add_chain(segment:Chain,new_id:str,structure:Structure)->None:
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
    chain_for_extraction=sys.argv[2]
    structure=read_in(test_path)
    output=extract_hetatm(structure[0][chain_for_extraction],inplace=True)
    add_chain(output[1],'Z',structure[0])
    write_out(structure)