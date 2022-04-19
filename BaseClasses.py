from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.Data import SCOPData
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import pandas as pd
import typing
from util import to_resid,read_in

allowed_level=typing.Literal['M','C','R']

class StructFeatureExtractor:
    '''
    '''
    def __init__(self)->None:
        '''
        '''
        self.operation_name=''
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
    
    def set_structure(self,struct:typing.Union[str,Structure])->None:
        '''
        '''
        if isinstance(struct,str):
            self.object=read_in(struct,struct)
        elif isinstance(struct,Structure):
            self.object=struct.copy()
        else:
            raise TypeError

        self.frame=pd.DataFrame()
        self._produce_dataframe()

    def _produce_dataframe(self)->None:
        '''
        '''
        raise NotImplementedError

    def transform(self)->None:
        raise NotImplementedError

    def extract_feature(self)->None:
        '''
        the structure can carry some feature in its`.xtra` property.
        use this method to put them into the frame.
        '''
        raise NotImplementedError
                
    def get_feature(self,feature_name:str)->pd.Series:
        return self.frame[feature_name].copy(deep=True)

class ModelFeatureExtractor(StructFeatureExtractor):
    '''
    '''
    def _produce_dataframe(self)->None:
        '''
        '''
        index_list_=[]
        for model in self.object:
            index_list_.append(model.id)
        index_=pd.Index(index_list_,name=['model'])
        self.frame=pd.DataFrame(index=index_)
    

class ChainFeatureExtractor(StructFeatureExtractor):  
    '''
    ''' 
    def _produce_dataframe(self)->None:
        index_list_=[[],[]]
        for model in self.object:
            for chain in model:
                index_list_[0].append(model.id)
                index_list_[1].append(chain.id)
        index_=pd.MultiIndex.from_arrays(index_list_,names=['model','chain'])
        self.frame=pd.DataFrame(index=index_)

class ResidueFeatureExtractor(StructFeatureExtractor):
    '''
    ''' 
    def _produce_dataframe(self)->None:
        '''
        '''
        index_list_=[[],[],[]]
        for model in self.object:
            for chain in model:
                for residue in chain:
                    index_list_[0].append(model.id)
                    index_list_[1].append(chain.id)
                    index_list_[2].append(residue.id)
        index_=pd.MultiIndex.from_arrays(index_list_,names=['model','chain','residue'])
        self.frame=pd.DataFrame(index=index_)
        return index_

    


class DistanceRelyFeatureExtractor(StructFeatureExtractor):
    pass



class BatchProcessor:
    pass