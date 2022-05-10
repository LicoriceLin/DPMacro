from Bio.PDB.StructureBuilder import StructureBuilder
# from Bio.Data import SCOPData
# from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
import pandas as pd
import numpy as np
from typing import Union,Literal,Tuple,Dict
from util import to_resid,read_in,extract_hetatm
from Data import amino3to1dict
# allowed_level=typing.Literal['M','C','R']

class StructFeatureExtractor:
    '''
    '''
    def __init__(self,operation_name:str='',canonical_only:bool=True)->None:
        '''
        '''
        self.operation_name=operation_name
        self._canonical_only=canonical_only
        self.object= Structure('')
        self.frame=pd.DataFrame()
        pass

    def __repr__(self):
        '''
        '''
        if hasattr(self,'object'):
            struct_id_=self.object.id
        else:
            struct_id_='null'
        if hasattr(self,'frame'):
            features_=', '.join(list(self.frame.columns))
        else:
            features_='null'
        return f'<Extractor : {self.operation_name} : {struct_id_} : {features_}>'
    
    def _set_structure(self,struct:Union[str,Structure])->None:
        '''
        '''
        if isinstance(struct,str):
            self.object=read_in(struct,struct)
        elif isinstance(struct,Structure):
            self.object=struct.copy()
        else:
            raise TypeError
        
    def _remove_hetatm(self)->None:
        '''
        '''
        for _chain in self.object.get_chains():
            extract_hetatm(_chain,inplace=True)

    def _init_dataframe(self)->None:
        '''
        '''
        raise NotImplementedError

    def _produce_feature(self)->None:
        '''
        '''
        pass

    def transform(self,struct:Union[str,Structure])->pd.DataFrame:
        '''
        '''
        self._set_structure(struct)
        if self._canonical_only:
            self._remove_hetatm()
        self._init_dataframe()
        self._produce_feature()
        return self.frame

    def _object_feature_to_frame(self)->None:
        '''
        the structure can carry features in its`.xtra` property,
        use this method to put them into the frame.
        '''
        raise NotImplementedError

    def _frame_feature_to_object(self)->None:
        '''
        reverse to `_object_feature_to_frame`
        '''
        raise NotImplementedError
                
    def get_feature(self)->pd.Series:
        '''
        the MultiIndex of pandas may be too tricky.
        I'm trying to wrap it into a more convenient manner  
        '''
        raise NotImplementedError


class ModelFeatureExtractor(StructFeatureExtractor):
    '''
    '''
    def _init_dataframe(self)->None:
        '''
        '''
        _index_list=[]
        for model in self.object:
            _index_list.append(model.id)
        _index=pd.Index(_index_list,name=['model'])
        self.frame=pd.DataFrame(index=_index)
    

class ChainFeatureExtractor(StructFeatureExtractor):  
    '''
    ''' 
    def _init_dataframe(self,from_residue:bool=False)->None:
        _index_list=[]
        for chain in self.object.get_chains():
            _index_list.append(chain.get_full_id()[1:])
        # _index_array=np.array(_index_list,dtype=object).transpose(1,0)
        # _index=pd.MultiIndex.from_arrays(_index_array,names=['model','chain'])
        self.frame=pd.DataFrame(_index_list,columns=['model','chain'])

        if from_residue:
            _index_list=[]
            for residue in self.object.get_residues():
                _index_list.append(residue.get_full_id()[1:])
            self.residue_frame=pd.DataFrame(_index_list,columns=['model','chain','residue'])

    
    def get_feature(self,feature,frame,chain)->pd.Series:
        '''
        '''
        return self.frame.loc[frame,chain][feature]


class ResidueFeatureExtractor(StructFeatureExtractor):
    '''
    ''' 
    def _init_dataframe(self)->None:
        '''
        '''
        _index_list=[]
        object_list=[]
        for residue in self.object.get_residues():
            _index_list.append(residue.get_full_id()[1:]+(residue.get_resname(),))
            object_list.append(residue)
        # _index_array=np.array(_index_list,dtype=object).transpose(1,0)
        # _index=pd.MultiIndex.from_arrays(_index_array,names=['model','chain','residue'])
        # self.frame=pd.DataFrame(index=_index)
        self.frame=pd.DataFrame(_index_list,columns=['model','chain','residue','resname'])
        self.frame['object']=object_list
    def _object_feature_to_frame(self) -> None:
        '''
        '''
        _init_flag=True
        for _residue in self.object.get_residues():
            _res_feats=_residue.xtra
            if _init_flag:
                _feature_dicts={key:[] for key in _res_feats.keys()}
                _init_flag=False
            for key,value in _res_feats.items():
                _feature_dicts[key].append(value)
        
        for key,value in _feature_dicts.items():
            self.frame[key]=value

            


class DistanceRelyFeatureExtractor(StructFeatureExtractor):
    '''
    '''
    def _produce_tree(self,entity:Entity)->None:
        pass

    def _distance_between_entity(self,
    entity1:Entity,entity2:Entity)->Tuple[tuple,tuple,float]:
        pass


class SeqFeatureExtractor:
    '''
    '''
    def __init__(self):
        pass

    def _set_seqs(self,seqs:str):
        self.seq=seqs
        pass

    def _init_dataframe(self):
        _seq_list=list(self.seq)
        _resname=[amino3to1dict.get(i,'XXX') for i in _seq_list]
        self.frame=pd.DataFrame(_seq_list,columns=['resname1'])
        self.frame['resname']=_resname
    def _produce_feature(self)->None:
        pass

    def transform(self,input:str)->pd.DataFrame:
        self._set_seqs(input)
        self._init_dataframe()
        self._produce_feature()
        return self.frame

class PretrainedFeatureExtractor:
    pass
class BatchProcessor:
    pass
class DirtyTmpfileFeatureExtractor:
    pass