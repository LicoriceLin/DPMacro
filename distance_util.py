'''
useful tools involving spacial calculation.
'''

from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
from Bio.PDB.Residue import Residue
from typing import List,Tuple,Union,Dict
# import typing
from Bio.PDB.Atom import Atom
from Bio.PDB.kdtrees import KDTree
import numpy as np
import pandas as pd
from util import integrated_residue_iterator,integrated_atom_iterator,allowed_residue_source
# import warnings
# from Bio import BiopythonParserWarning




#true functions


def produce_tree(object:allowed_residue_source)->Tuple[pd.DataFrame,KDTree]:
    '''
    yield the basic instance for distance calculation.
    output[0]: a dataframe:
        'model','chain','residue','atom','resname':atom information;
        'x','y','z': coordinates;
        'object': Atom instance. extact the same instance as Atoms in input.
    output[1]: a kdtree object implemented in biopython 1.79+.
        use `from Bio.PDB.kdtrees import KDTree ; help(KDTree) ` to see the detail.
    '''
    _atom_index_list=[]
    _atom_coord_list=[]
    _atom_resname_list=[]
    _atom_list=[]
    for atom in integrated_atom_iterator(object):
        _full_id=atom.get_full_id()[1:]
        _atom_index_list.append([idtuple2str(i) for i in _full_id])
        _atom_coord_list.append(atom.coord)
        _atom_resname_list.append([atom.get_parent().get_resname()])
        _atom_list.append(atom)
    
    atom_index=np.array(_atom_index_list,dtype=str)
    atom_coord=np.array(_atom_coord_list,dtype=np.float64)
    atom_resname=np.array(_atom_resname_list,dtype=str)
    tree=KDTree(atom_coord,10)

    frame=pd.DataFrame(np.concatenate([atom_index,atom_resname],axis=1),columns=['model','chain','residue','atom','resname'])
    coord_frame=pd.DataFrame(atom_coord,columns=list('xyz'))
    atom_frame=pd.DataFrame(_atom_list,columns=['object'])
    frame=frame.join(coord_frame,how='right')
    frame=frame.join(atom_frame,how='right')
    return frame,tree

def distance_between_entity(
        entity1:Union[Entity,List[Entity]],entity2:Union[Entity,List[Entity]])->Tuple[Atom,Atom,float]:
    '''
    the closest distance between 2 atom sets.
    return[0] / return[1] : the closest Atom instance in entity1 & entity2.
    return[2] : the distance.
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

def atom_within_threshold(
        entity1:Union[Entity,List[Entity]],entity2:Union[Entity,List[Entity]],threshold:float)->Dict[Atom,float]:
    '''
    search entity1, 
    return {Atom instance within the threshold of entity2:the distance}
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

def residue_within_threshold(entity1:allowed_residue_source,entity2:allowed_residue_source,threshold)->List[Residue]:
    '''
    return a non-redundant list of Residue in entity1 within distance threshold of entity2.
    Note:only search atom in entity1, get atom within threshold,
        and then get their parent in a deduplicated list.
    '''
    return list(set(integrated_residue_iterator(
                            atom_within_threshold(entity1,entity2,threshold).keys())))

def residue_distance_matrix(object:allowed_residue_source)->pd.Series:
    '''
    return a residue-wise closest distance matrix from the input Atom set.
    use output[Residue1,Residue2] to get the distance
    note: if only part of a residue's atom is included in the input, clalulation will consider only these atoms.   

    '''
    frame=produce_tree(object)[0][['object','x','y','z']]
    frame['residue']=frame['object'].apply(Atom.get_parent)
    cross_frame=get_distance(frame,frame)
    distances = cross_frame.groupby(['residue_1','residue_2'])['distance'].min()
    return distances

def atom_distance_matrix(object:allowed_residue_source)->pd.Series:
    '''
    return a atom-wise closest distance matrix from the input Atom set.
    use output[Atom1,Atom2] to get the distance between Atom1&2

    '''
    frame=produce_tree(object)[0][['object','x','y','z']]
    frame.columns=['atom','x','y','z']
    # frame['residue']=frame['object'].apply(Atom.get_parent)
    cross_frame=get_distance(frame,frame)
    distances = cross_frame['distance']
    distances.index=pd.Index(cross_frame[['atom_1','atom_2']])
    return distances

#small tools
def idtuple2str(idtuple:Union[Tuple,str,int])->str:
    '''
    (deprecated)
    turn "standard triplet biopython residue code" to a single string linked by "|"
    '''
    if isinstance(idtuple,tuple):
        idlist=[str(i) if i != ' ' else '_'  for i in idtuple  ]
        return '|'.join(idlist)
    else:
        return str(idtuple) if idtuple != ' ' else '_'
    
def idstr2tuple(idstr:str)->Tuple:
    '''
    (deprecated)
    turn  a single string linked by "|" to "standard triplet biopython residue code"
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
    (deprecated)
    used in produce_tree's dataframe.
    extract the Atom's full_id from this dataframe.
    '''
    assert set(['model','chain','residue','atom'])<set(frame.columns),'frame does not contain adequate columns'
    _id=frame.loc[index][['model','chain','residue','atom']]
    id=tuple(_id.map(idstr2tuple))
    return id

def atom_object_from_frame(frame:pd.DataFrame,index:int)->Atom:
    '''
    (actually in use )
    used in produce_tree's dataframe.
    extract the Atom object from dataframe.
    '''
    assert set(['object'])<set(frame.columns),'frame does not contain adequate columns'
    object=frame.loc[index]['object']
    return object

def residue_id_from_frame(frame:pd.DataFrame,index:int)->Tuple[int,str,tuple,tuple]:
    '''
    (deprecated)
    used in produce_tree's dataframe.
    extract the Atom.get_parent()'s full_id, where the Atom's index in the dataframe ( and the tree) is given.
    '''
    assert set(['model','chain','residue'])<set(frame.columns),'frame does not contain adequate columns'
    _id=frame.loc[index][['model','chain','residue']]
    id=tuple(_id.map(idstr2tuple))
    return id

def _euclid_dis(x1:float,y1:float,z1:float,x2:float,y2:float,z2:float)->float:
    '''
    calculate the spacial distance.
    '''
    return  ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5

def get_distance(res1:pd.DataFrame,res2:pd.DataFrame)->pd.DataFrame:
    '''
    input:2 dataframes from produce_tree
    return a distance matrix, containing all distance between atoms in res1 and res2.
    use return[(return['object_1']==Atom1) & (return['object_2']==Atom2)]['distance'] to get distance between Atom1 in frame1 and Atom2 in frame2
    '''
    dis=pd.merge(res1,res2,how='cross',suffixes=('_1','_2'))
    dis['distance']=((dis['x_1']-dis['x_2'])**2
        +(dis['y_1']-dis['y_2'])**2
        +(dis['z_1']-dis['z_2'])**2)**0.5
    return dis


#deprecate functions
def unintegrate_produce_tree(entity:Union[Entity,Atom,List[Entity],List[Atom]])->Tuple[pd.DataFrame,KDTree]:
    '''
    (deprecated)
    olde version of produce_tree
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

def frame_distance_between_entity(
        entity1:Union[Entity,List[Entity]],entity2:Union[Entity,List[Entity]])->Tuple[Tuple,Tuple,float]:
    '''
    (deprecated)
    slow version of `distance_between_entity`
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

def tree_distance_matrix(struct:Structure)->pd.DataFrame:
    '''
    (deprecated)
    slow version of `distance_matrix`
    '''
    # if len(object.child_list)>1:
    #     warnings.warn("multi-frame object detected."
    #          "only frame 0 will be processed.",
    #          BiopythonParserWarning,)
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


#need update
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
