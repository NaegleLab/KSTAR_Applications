import pickle
import os
import pandas as pd


"""

This is where you can set the GLOBAL Path directory used in examples to point to the supplementary data for KSTAR

"""

SUPPLEMENTS_DIR = '../../data/KSTAR_NETWORKS/KSTAR_July2021/SupplementaryData/'

KINASE_MAP =  pd.read_csv(SUPPLEMENTS_DIR+'Map/globalKinaseMap.csv', index_col = 0)
