from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.Data import SCOPData
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import pandas as pd
import typing
from util import to_resid,read_in,extract_hetatm

allowed_level=typing.Literal['M','C','R']

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
    
    def _set_structure(self,struct:typing.Union[str,Structure])->None:
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

    def transform(self,struct:typing.Union[str,Structure])->pd.DataFrame:
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
    def _init_dataframe(self)->None:
        _index_list=[[],[]]
        for model in self.object:
            for chain in model:
                _index_list[0].append(model.id)
                _index_list[1].append(chain.id)
        _index=pd.MultiIndex.from_arrays(_index_list,names=['model','chain'])
        self.frame=pd.DataFrame(index=_index)

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
        _index_list=[[],[],[]]
        for model in self.object:
            for chain in model:
                for residue in chain:
                    _index_list[0].append(model.id)
                    _index_list[1].append(chain.id)
                    _index_list[2].append(residue.id)
        _index=pd.MultiIndex.from_arrays(_index_list,names=['model','chain','residue'])
        self.frame=pd.DataFrame(index=_index)

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


class DistanceRelyFeatureExtractor:
    def __init__(self)->None:
        pass
    def _produce_tree(self)->None:
        pass

class PretrainedFeatureExtractor:
    pass

class DirtyTmpfileFeatureExtractor:
    pass

class SeqFeatureExtractor:
    pass

class BatchProcessor:
    pass