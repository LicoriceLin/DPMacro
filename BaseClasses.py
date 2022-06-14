'''
BaseClasses
maintains interfaces and shared method.
        '''



from Bio.PDB.StructureBuilder import StructureBuilder
# from Bio.Data import SCOPData
# from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
import pandas as pd
# import numpy as np
from typing import Union,Tuple,Dict,Iterable
from util import read_in,extract_hetatm,sequence_from_frame
from util import _list_feature_into_residue,_list_feature_into_frame,integrated_residue_iterator
from util import allowed_residue_source,allowd_scheme,RESIDUE_HOLDER
from Data import amino3to1dict

class StructFeatureExtractor:
    '''
    parent of Model/Chain/ResidueFeatureExtractor
    but the first 2 haven't been implemented yet.  

    public property:
    self.object
        Bio.Structure instance.
        some quick guide:
            1. In Bio. pdb information is maintained hierarchical in `Structure`(like Universe in mdanalyse), `Model`(frame), Chain, Residue and Atom. they are all subclass of `Entity` except for Atom.
            2. Entity is Iterable. e.g. 
                Structure[0] returns the first frame in the pdb file; 
                Model['A'] returns chain A in the frame; 
                Chain[10] or Chain[(' ',107,'A')](use the latter to get residue with insertion-code) returns a Residue
                Residue['CB'] returns a Atom.
            3. Both Entity and Atom are Hashable.
            4.I use a dict, Entity/Atom.xtra to save features in the structure in the object.
            See More in the Biopython Tutorial.
        
    self.frame
        pd.DataFrame object.
        See more in subclasses of StructFeatureExtractor

    self.operation_name:
        a str to display ( and maybe hash, in the future) this extractor instance.

    method:
        see docstring of each method.

    '''
    def __init__(self,operation_name:str='',canonical_only:bool=True,if_reduce:bool=False)->None:
        '''
        set the name of your Extractor
        decide the parameters for transform. 
        (e.g. if solvent and ligand will be maintained in the structure) 
        '''
        self.operation_name=operation_name
        self._canonical_only=canonical_only
        self._if_reduce=if_reduce
        self.object= Structure('')
        self.frame=pd.DataFrame()
        pass

    def __repr__(self):
        '''
        show Extractor instance in following manner:
        <Extractor : {self.operation_name} : {struct_id_} : ', '.join(list(self.frame.columns))>
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
        see `transform(self)`
        '''
        if isinstance(struct,str):
            self.object=read_in(struct,struct)
        elif isinstance(struct,Structure):
            self.object=struct.copy()
        else:
            raise TypeError
        if self._canonical_only:
            self._remove_hetatm()
        
    def _remove_hetatm(self)->None:
        '''
        delete all the ligand and solvent in self.object. 
        warning: untest in ncaa systems.

        '''
        for _chain in self.object.get_chains():
            extract_hetatm(_chain,inplace=True)

    def _init_dataframe(self)->None:
        '''
         (they are just interface to be implemented in subclass.)
        see their expected behavior in `transform(self)`
        '''
        raise NotImplementedError

    def _produce_feature(self)->None:
        '''
        (interface to be implemented in subclass.)
        see their expected behavior in `transform(self)`
        '''
        pass

    def transform(self,struct:Union[str,Structure],**kwargs)->Union[pd.DataFrame,Dict,Iterable]:
        '''
        1.set the input structure( both Bio.Structure object and path to pdb file is allowed) to self.object
        2.do some simple pre-process( e.g. remove hetatm, split Fv domain)
        3.init self.frame as a dataframe containing all entity with feature( e.g. a frame in which each row points to a residue )
        4.produce the feature by private method `_produce_feature`, note all the parameter for this method should be defined during `__init__` of the instance.
        5.return self.frame. contains entity and their features.
        '''
        self._set_structure(struct)
        self._init_dataframe()
        self._produce_feature(**kwargs)
        if self._if_reduce:
            self._reduce()
            return self.reduced_frame
        else:
            return self.frame

    def _object_feature_to_frame(self)->None:
        '''
        (interface to be implemented in subclass.)
        the structure can carry features in its`.xtra` property,
        use this method to put them into the frame.
        refresh self.frame by pull features from self.object.
        note: supposed to work in an `append new and overwrite old` manner 
        '''
        raise NotImplementedError

    def _frame_feature_to_object(self)->None:
        '''
        (interface to be implemented in subclass.)
        reverse to `_object_feature_to_frame`
        refresh self.object by pull features from self.frame.
        note: supposed to work in an `append new and overwrite old` manner 
        '''
        raise NotImplementedError
                
    def get_feature(self)->pd.Series:
        '''
        may be deprecated soon
        the MultiIndex of pandas may be too tricky.
        I'm trying to wrap it into a more convenient manner  
        '''
        raise NotImplementedError

    def _reduce(self):
        '''
        (interface for rewrite in its subclass.)
        reduce Residue property to Chain-level or Model-level property
        expected behavior:
        deal with the self.frame.groupby(...) object.
        return a chain or model level feature frame.
        '''
        self.reduced_frame:Iterable=pd.DataFrame()

class ResidueFeatureExtractor(StructFeatureExtractor):
    '''
    subclass of StructFeatureExtractor.
    used to process Residue-level features, 
    and reduce them to Chain/Model level
    ''' 
    def _init_dataframe(self)->None:
        '''
        note the format in ResidueFeatureExtractor.frame

        a fresh initiated self.frame contains 5 columns:
            model: model id for each residue;
            chain: chain id for each residue;
            residue: residue `full id`(a trimer tuple, see biopython tutorial for detail) for each residue;
            resname: residue's full name (three letter, all capital)
            object: the Residue object itself. 
                !!it is the SAME instance as Residue in self.object.
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
        see parent.
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

    def _produce_sequence(self,map:Dict[str,str]=amino3to1dict)->None:
        # def list3_to_str1(list3:Iterable[str])->str:
        #     return ''.join([map.get(i,'X') for i in list3])
        # _seq_series=self.frame.groupby('chain')['resname'].apply(list3_to_str1)
        # self.sequences=dict(_seq_series)
        '''
        initiate a {chain:fastaseq} dict in self.sequences by running util.sequence_from_frame
        '''
        self.sequences=sequence_from_frame(self.frame,map)



#unimplemented
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

class SeqFeatureExtractor:
    '''
    under early construction.
    share similar behavior with StructFeatureExtractor;
    but use a sequence str as input.
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


#may be deprecated
# class DistanceRelyFeatureExtractor(StructFeatureExtractor):
#     '''
#     may be deprecated
#     '''
#     def _produce_tree(self,entity:Entity)->None:
#         pass

#     def _distance_between_entity(self,
#     entity1:Entity,entity2:Entity)->Tuple[tuple,tuple,float]:
#         pass

# class PretrainedFeatureExtractor:
#     '''
#     under development
#     '''
#     pass

# class BatchProcessor:
#     '''
#     may be deprecated
#     '''
#     pass

# class DirtyTmpfileFeatureExtractor:
#     '''
#     may be deprecated
#     '''
#     pass