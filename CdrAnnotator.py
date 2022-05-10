from pandas import DataFrame
import distance_util as du
import pandas as pd
import numpy as np
from BaseClasses import ResidueFeatureExtractor,ChainFeatureExtractor,SeqFeatureExtractor
from Bio.PDB.Entity import Entity
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Data import *
from Bio import BiopythonParserWarning
from typing import List,Tuple,Union,Dict
import warnings