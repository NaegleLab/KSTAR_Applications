import pickle
import os
import pandas as pd


"""

This is where you can set the GLOBAL Path directory used in examples to point to the supplementary data for KSTAR

"""

SUPPLEMENTS_DIR = '/Users/kmn4mj/Box Sync/Manuscripts/2021/KSTAR/'

KINASE_MAP =  pd.read_csv(SUPPLEMENTS_DIR+'Supplements/SupplementaryData/Map/globalKinaseMap.csv', index_col = 0)
