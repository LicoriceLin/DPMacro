__version__ = '0.3.0dev0'
__author__ = "Zhenfeng Deng"
import warnings

try:
    from .util import (integrated_residue_iterator,read_in,write_out,impute_beta,
                    split_frame,extract_hetatm,sequence_from_frame)
    from .distance_util import (
        atom_distance_matrix, residue_distance_matrix,
        atom_within_threshold, residue_within_threshold,
        distance_between_entity
        
    )
except:
    warnings.warn(('DPMacro fundament'
                   ' does not init properly!'))

try:
    from .FvStructureProcesser import FvStructureProcesser,FvSpliter
    from .CdrAnnotator import hmt_FvProcessor
except:
    warnings.warn(('Antibody module'
                   ' does not init properly!'))
    
try:
    from .AFill import structure_aligner,reslist_to_guicode
    from .PockerExtractor import PocketExtractor
except:
    warnings.warn(('Pocket module'
                   ' does not init properly!'))
    
try:
    from .TAPFeature import TAPExtractor
    from .PpsFeatureExtractor import PPSExtractor
except:
    warnings.warn(('Feature module'
                   ' does not init properly!'))


# __all__ = ["FvStructureProcesser","FvSpliter"]