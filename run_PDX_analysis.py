#Import preamble of kstar and other necessary functions
import pandas as pd
import os
import pickle

from kstar import config, helpers, calculate


import conf #applications config for supplementary directory

data_dir = conf.SUPPLEMENTS_DIR+'BreastCancer/PDX_Huang2017/'

odir = './data/predictions/'
expName = 'PDX'

experiment = pd.read_csv(f"{data_dir}/MAPPED/{expName}_mapped.tsv", sep='\t')

networks = {}
networks['Y'] = pickle.load(open(config.NETWORK_Y_PICKLE, "rb" ) )


if not os.path.exists(f"{odir}/{expName}"): 
    os.mkdir(f"{odir}/{expName}")
activity_log = helpers.get_logger(f"activity_{expName}", f"{odir}/{expName}/activity_{expName}.log")

data_columns = None # this means default to the column headers marked with data:
agg = 'mean' # if a non-NaN value appears at all, use it
threshold = 0.0 #require they have increased from basal

#Run analysis, normalization, and Mann Whitney. 
#Save the SLIM versions at each step, in case of failure/loss of kernel, etc.

kinact_dict = calculate.run_kstar_analysis(experiment, activity_log, networks, phospho_types = ['Y'], data_columns = data_columns, agg =agg, threshold = threshold,  greater = True)
calculate.save_kstar_slim(kinact_dict, expName, f"{odir}/{expName}")


num_random_experiments=150
target_alpha=0.05
calculate.normalize_analysis(kinact_dict, activity_log, num_random_experiments, target_alpha)
calculate.save_kstar_slim(kinact_dict, expName, f"{odir}/{expName}")


kinact_dict['Y'].calculate_Mann_Whitney_activities_sig(activity_log, number_sig_trials = 100)
calculate.save_kstar_slim(kinact_dict, expName, f"{odir}/{expName}")