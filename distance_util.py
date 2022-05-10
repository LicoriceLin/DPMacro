from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
from typing import List,Tuple,Union,Dict
import typing
from Bio.PDB.Atom import Atom
from Bio.PDB.kdtrees import KDTree
import numpy as np
import pandas as pd
from util import idtuple2str,idstr2tuple,atom_object_from_frame
import warnings
from Bio import BiopythonParserWarning


def produce_tree(entity:Union[Entity,Atom,List[Entity],List[Atom]])->Tuple[pd.DataFrame,KDTree]:
    '''
    '''
    _atom_index_list=[]
    _atom_coord_list=[]
    _atom_resname_list=[]
    _atom_list=[]
    if isinstance(entity,list):
        for _sub_entity in entity:
            if not isinstance(_sub_entity,Atom):
                for _atom in _sub_entity.get_atoms():
                    _full_id=_atom.get_full_id()[1:]
                    _atom_index_list.append([idtuple2str(i) for i in _full_id])
                    _atom_coord_list.append(_atom.coord)
                    _atom_resname_list.append([_atom.get_parent().get_resname()])
                    _atom_list.append(_atom)
            else:
                _atom=_sub_entity
                _full_id=_atom.get_full_id()[1:]
                _atom_index_list.append([idtuple2str(i) for i in _full_id])
                _atom_coord_list.append(_atom.coord)
                _atom_resname_list.append([_atom.get_parent().get_resname()])
                _atom_list.append(_atom)
    elif isinstance(entity,Entity):
        for _atom in entity.get_atoms():
            _full_id=_atom.get_full_id()[1:]
            _atom_index_list.append([idtuple2str(i) for i in _full_id])
            _atom_coord_list.append(_atom.coord)
            _atom_resname_list.append([_atom.get_parent().get_resname()])
            _atom_list.append(_atom)
    elif isinstance(entity,Atom):
        _atom=entity
        _full_id=_atom.get_full_id()[1:]
        _atom_index_list.append([idtuple2str(i) for i in _full_id])
        _atom_coord_list.append(_atom.coord)
        _atom_resname_list.append([_atom.get_parent().get_resname()])
        _atom_list.append(_atom)
    else:
        raise TypeError
    
    atom_index=np.array(_atom_index_list,dtype=str)
    atom_coord=np.array(_atom_coord_list,dtype=np.float64)
    atom_resname=np.array(_atom_resname_list,dtype=str)
    tree=KDTree(atom_coord,10)

    # full_array=np.concatenate([_atom_index_list,_atom_resname_list,_atom_coord_list],axis=1).T
    # _index=pd.MultiIndex.from_arrays(atom_index.transpose(1,0),names=['model','chain','residue','atom'])
    # frame=pd.DataFrame(atom_coord,index=_index,columns=['x','y','z'])
    frame=pd.DataFrame(np.concatenate([atom_index,atom_resname],axis=1),columns=['model','chain','residue','atom','resname'])
    coord_frame=pd.DataFrame(atom_coord,columns=list('xyz'))
    atom_frame=pd.DataFrame(_atom_list,columns=['object'])
    frame=frame.join(coord_frame,how='right')
    frame=frame.join(atom_frame,how='right')
    # return atom_index,atom_coord,atom_resname,frame.join(coord_frame,how='right'),tree
    return frame,tree

def distance_between_entity(
        entity1:Union[Entity,List[Entity]],entity2:Union[Entity,List[Entity]])->Tuple[Atom,Atom,float]:
    '''
    '''
    frame1,tree1=produce_tree(entity1)
    frame2,tree2=produce_tree(entity2)
    distance=9999.9
    entity1_index='_'
    entity2_index='_'
    for i in frame1[['x','y','z']].iterrows():
        _entity1_index=i[0]
        entity1_cord=np.array(i[1],dtype=np.float64)
        search=tree2.search(entity1_cord,1000)
        _min=np.argmin([x.radius for x in search])
        _entity2_index=search[_min].index
        _distance=search[_min].radius
        # print(f'{distance:.2f}',end='\t')
        if _distance<distance:
            distance=_distance
            entity1_index=_entity1_index
            entity2_index=_entity2_index
    # _id1=frame1.loc[entity1_index][['model','chain','residue','atom']]
    # _id2=frame2.loc[entity2_index][['model','chain','residue','atom']]
    # id1=tuple(_id1.map(idstr2tuple))
    # id2=tuple(_id2.map(idstr2tuple))
    id1=atom_object_from_frame(frame1,entity1_index)
    id2=atom_object_from_frame(frame2,entity2_index)

    return id1,id2,distance

def _euclid_dis(x1:float,y1:float,z1:float,x2:float,y2:float,z2:float):
    '''
    '''
    return  ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5

def frame_distance_between_entity(
        entity1:Union[Entity,List[Entity]],entity2:Union[Entity,List[Entity]])->Tuple[Tuple,Tuple,float]:
    '''
    much slower
    '''
    frame1,tree1=produce_tree(entity1)
    frame2,tree2=produce_tree(entity2)
    crossframe=pd.merge(left=frame1,right=frame2,how='cross',suffixes=('_1','_2'))
    crossframe['distance']=np.vectorize(_euclid_dis)(crossframe['x_1'],crossframe['y_1'],crossframe['z_1'],
                                            crossframe['x_2'],crossframe['y_2'],crossframe['z_2'])
    near_idx=crossframe['distance'].idxmin()
    distance=crossframe.loc[near_idx]['distance']
    _id1=crossframe.loc[near_idx][['model_1','chain_1','residue_1','atom_1']]
    _id2=crossframe.loc[near_idx][['model_2','chain_2','residue_2','atom_2']]
    id1=tuple(_id1.map(idstr2tuple))
    id2=tuple(_id2.map(idstr2tuple))
    return id1,id2,distance

def atom_within_threshold(
        entity1:Union[Entity,List[Entity]],entity2:Union[Entity,List[Entity]],threshold:float)->Dict[Tuple,float]:
    '''
    '''
    frame1,tree1=produce_tree(entity1)
    frame2,tree2=produce_tree(entity2)
    _outputdict={}
    for i in frame2[['x','y','z']].iterrows():
        entity2_cord=np.array(i[1],dtype=np.float64)
        for i in tree1.search(entity2_cord,threshold):
            _outputdict[i.index]=min(i.radius,_outputdict.get(i.index,9999.9))
    outputdict={atom_object_from_frame(frame1,key):value for key,value in _outputdict.items()}
    return outputdict

def residue_distance_matrix(struct:Structure)->pd.DataFrame:
    if len(object.child_list)>1:
        warnings.warn("multi-frame object detected."
             "only frame 0 will be processed.",
             BiopythonParserWarning,)
    index=[]
    dis_matrix=[]
    for residue in struct[0].get_residues():
        index.append(residue)
        dis_line=[]
        for residue_ in struct.get_residues():
            dis_line.append(distance_between_entity(residue,residue_)[2])
        dis_matrix.append(dis_matrix)
    dis_frame=pd.DataFrame(dis_matrix,index=index,columns=index)
    return dis_frame

def get_distance(res1:pd.DataFrame,res2:pd.DataFrame)->float:
    dis=pd.merge(res1,res2,how='cross',suffixes=('_1','_2'))
    dis['distance']=((dis['x_1']-dis['x_2'])**2
        +(dis['y_1']-dis['y_2'])**2
        +(dis['z_1']-dis['z_2'])**2)**0.5
    return dis

def frame_residue_distance_matrix(struct:Structure)->pd.Series:
    frame=produce_tree(struct)[0][['object','x','y','z']]
    frame['residue']=frame['object'].apply(Atom.get_parent)
    cross_frame=get_distance(frame,frame)
    distances = cross_frame.groupby(['residue_1','residue_2'])['distance'].min()
    return distances


    

if __name__ == '__main__':
    import sys
    from util import read_in
    testfile=sys.argv[1]
    testcase=read_in(testfile)
    chain1=sys.argv[2]
    chain2=sys.argv[3]
    chain3=sys.argv[4]
    print(distance_between_entity(testcase[0][chain1],[testcase[0][chain2],testcase[0][chain3]]))
    print(atom_within_threshold(testcase[0][chain1],[testcase[0][chain2],testcase[0][chain3]],5))
