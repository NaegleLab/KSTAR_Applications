import pickle
import os
import pandas as pd


"""

This is where you can set the GLOBAL Path directory used in examples to point to the supplementary data for KSTAR, 
which can be found in the KSTAR Figshare.

"""
#only change this line to point to FigShare data
SUPPLEMENTS_DIR = './FigShare_Rev3_March30_2022/'

KINASE_MAP =  pd.read_csv(SUPPLEMENTS_DIR+'Map/globalKinaseMap.csv', index_col = 0)
