#benchmarking Kinase Library Python tool 
###https://github.com/TheKinaseLibrary/kinase-library
#using benchmarking data from KSTAR paper
#10.1038/s41467-022-32017-5
#re-preprocessed and re-mapped March 2026 by Sam Crowl
#https://myuva.sharepoint.com/:f:/r/sites/NaegleLab/Shared%20Documents/General/Projects/PhosphoproteomicsDatabase/KSTAR_Benchmarking/KinaseActivity_Benchmarks?csf=1&web=1&e=LCg6b1

#%%IMPORT PACKAGES
import os
from importlib import reload

import pandas as pd
import numpy as np

import kinase_library as kl
from kinase_library.modules.enrichment import plot_volcano
from kstar import config

import matplotlib.pyplot as plt
import seaborn as sns

import pickle

#%%FUNCTIONS FOR ANALYSIS
#import Sam's benchmarking functions
import benchmarking_functions

def convert_pep_to_central_format(peptide, pad=(7,7)):
    '''Convert peptide sequence to central format for KL enrichment analysis. Will infer phosphorylated amino acid from peptide sequence, lowercase letter indicates phosphorylated amino acid.
    Parameters:
    peptide: str
        peptide sequence in format 'AA(pAA)AA' where (pAA) indicates phosphorylated amino acid
    pad: tuple of ints
        number of amino acids to include on either side of the central phosphorylated amino acid (default is (7,7))
    Returns:
    central_format_peptide: str
        peptide sequence in format 'AAAAA[pAA]AAAAA' where [pAA] indicates phosphorylated amino acid and there are pad[0] amino acids on the left and pad[1] amino acids on the right. 
    '''
    #find index of phosphorylated amino acid
    pAA_options = ['s', 't', 'y']
    
    #iterate through possible phosphorylated amino acids
    #find index of phosphorylated amino acid and determine which amino acid is phosphorylated
    #assuming single amino acid is phosphorylated in peptide sequence
    for aa in pAA_options:
        if aa in peptide:
            p_index = peptide.index(aa)
            pAA = aa
            break
    
    #extract amino acids on either side of phosphorylated amino acid
    left_aa = peptide[:p_index]
    right_aa = peptide[p_index+len(f'{pAA}'):]
    
    #pad with _ if there are not enough amino acids on either side
    left_pad = '_' * max(0, pad[0] - len(left_aa))
    right_pad = '_' * max(0, pad[1] - len(right_aa))
    
    #construct central format peptide
    central_format_peptide = f'{left_pad}{left_aa}{pAA}{right_aa}{right_pad}'
    
    return central_format_peptide

def run_KL_diff_phos_enrichment_analysis(experiment, kin_type, 
                                         threshold=1, agg='mean', 
                                         reformat_peptides=True, 
                                         quant_type='fold_change',
                                         enrichment_type='both',
                                         rename_kinases=False, kinaseNames=None):
    '''Run KL differential phosphorylation enrichment analysis by providing df of phosphorylation data with peptide sequences and quantitative values. For single column of quantitative data.
    PARAMETERS:
    experiment: dataframe with phosphoproteomic data
        must contain columns for KSTAR_ACCESSION, KSTAR_SITE, KSTAR_PEPTIDE, and at least one data column with quantitative values (e.g. fold change)
    kin_type: str
        type of kinase to analyze ('tyrosine' or 'ser_thr')
    threshold: numeric
        should be same quant type as data columns
    agg: str
        aggregation method for duplicate sites ('mean', 'sum', False) (default: 'mean')
    reformat_peptides: bool
        whether to reformat peptides to central format for KL enrichment analysis (default: True)
    enrichment_type: str
        type of enrichment analysis to run ('both', 'enriched', 'depleted') (default: 'both')
    quant_type: str
        type of quantitative data in data columns ('fold_change', 'log2_fold_change', 'other') (default: 'fold_change')
    rename_kinases: bool
        whether to rename kinases in KL results to standard kinase names using kinaseNames resource file (default: False)   
    kinaseNames: dataframe with kinase names for mapping KL kinase names -> std names, must contain columns 'KL_name' and 'std_name' (required if rename_kinases=True)
    RETURNS:
    klRes: dataframe with KL enrichment results for data
    '''
    #REFORMAT EXPERIMENT
    data_col = [col for col in experiment.columns if col.startswith('data:')][0]
    #aggregate
    if agg:
        experiment_agg = experiment.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[data_col].agg(agg).reset_index()
        #add KSTAR_PEPTIDE column back to experiment
        seq_info = experiment[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE']]
        experiment = experiment_agg.merge(seq_info, on=['KSTAR_ACCESSION', 'KSTAR_SITE'], how='left')
    
    #reformat peptides to central format if reformat_peptides=True
    if reformat_peptides:
        experiment['KSTAR_PEPTIDE'] = experiment['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x))
        #capitalize peptides
        experiment['KSTAR_PEPTIDE'] = experiment['KSTAR_PEPTIDE'].str.upper()   

    #convert to log2 fold change if needed
    if quant_type == 'fold_change':
        #add a small constant to data to avoid log of 0
        experiment[data_col] = experiment[data_col].replace(0, 1e-10)
        experiment[data_col] = np.log2(experiment[data_col])
        lfc_thresh = np.log2(threshold)
    elif quant_type == 'log2_fold_change':
        lfc_thresh = threshold
    
    #CREATE DiffPhosData OBJECT
    #drop nans in data_col
    experiment = experiment.dropna(subset=[data_col])
    DiffPhosData = kl.DiffPhosData(dp_data=experiment, lfc_col=data_col, seq_col='KSTAR_PEPTIDE', 
                                   lfc_thresh=lfc_thresh)

    #RUN KL ENRICHMENT ANALYSIS
    #determine kl_thresh based on kin_type
    if kin_type == 'tyrosine':
        kl_thresh = 8
    elif kin_type == 'ser_thr':
        kl_thresh = 15
    #run analysis
    klRes = DiffPhosData.kinase_enrichment(kin_type=kin_type, kl_method='percentile_rank', 
                                   kl_thresh=kl_thresh, enrichment_type=enrichment_type, non_canonical=True)
    
    #REFORMAT RESULTS
    #add column for data column
    klRes = klRes.combined_enrichment_results
    klRes['data_column'] = data_col
    #add column for enrichment type
    #use Up/Down/Both - stay consistent with benchmarking meta dataframe
    enrichment_type_dict = {'enriched': 'Up', 'depleted': 'Down', 'both': 'Both'}
    enrichment_type = enrichment_type_dict.get(enrichment_type, enrichment_type)
    klRes['enrichment_type'] = enrichment_type
    #move index to column and name 'Kinase'
    klRes = klRes.reset_index().rename(columns={'index': 'kinase'})
    #rename kinases in klRes to standard names using kinaseNames resource file if rename_kinases=True
    if rename_kinases:
        if kinaseNames is None:
            raise ValueError('kinaseNames dataframe must be provided if rename_kinases=True')
        #make name mapping dict
        name_map = dict(zip(kinaseNames['KLETcst_py'], kinaseNames['ReName']))
        #store original names
        klRes['kinase_KLname'] = klRes['kinase'].copy()
        #map kinase names in klRes to standard names using name_map dict
        klRes['kinase'] = klRes['kinase_KLname'].map(name_map)
    
    return klRes

def run_KL_diff_phos_enrichment_analysis_batch(experiment, kin_type, 
                                         threshold=1, agg='mean', 
                                         reformat_peptides=True, 
                                         quant_type='fold_change',
                                         enrichment_type='both',
                                         meta=None, meta_data_col_name='benchmarking_column', meta_pertub_direction_col_name='perturbation_direction',
                                         rename_kinases=False, kinaseNames=None):
    '''Run KL differential phosphorylation enrichment analysis by providing df of phosphorylation data with peptide sequences and quantitative values. For single column of quantitative data.
    PARAMETERS:
    experiment: dataframe with phosphoproteomic data
        must contain columns for KSTAR_ACCESSION, KSTAR_SITE, KSTAR_PEPTIDE, and at least one data column with quantitative values (e.g. fold change)
    kin_type: str
        type of kinase to analyze ('tyrosine' or 'ser_thr')
    threshold: numeric
        should be same quant type as data columns
    agg: str
        aggregation method for duplicate sites ('mean', 'sum', False) (default: 'mean')
    reformat_peptides: bool
        whether to reformat peptides to central format for KL enrichment analysis (default: True)
    enrichment_type: str
        type of enrichment analysis to run ('both', 'enriched', 'depleted', None) (default: 'both')
    meta: dataframe with metadata for benchmarking dataset, must contain columns for dataset name and perturbation direction. Will override enrichment_type parameter to determine direction of enrichment analysis for each dataset in batch.
    meta_data_col_name: str
        column name in meta dataframe that contains column names in benchmarking dataset (should match the dataset names in the experiment dataframe) (default: 'benchmarking_column')
    meta_pertub_direction_col_name: str
        column name in meta dataframe that contains perturbation direction information (default: 'perturbation_direction')
    quant_type: str
        type of quantitative data in data columns ('fold_change', 'log2_fold_change', 'other') (default: 'fold_change')
    rename_kinases: bool
        whether to rename kinases in KL results to standard kinase names using kinaseNames resource file (default: False)   
    kinaseNames: dataframe with kinase names for mapping KL kinase names -> std names, must contain columns 'KL_name' and 'std_name' (required if rename_kinases=True)
    RETURNS:
    klRes: dataframe with KL enrichment results for data
    '''
   
    #RUN KL ENRICHMENT ANALYSIS
    #determine enrichment type for column based on meta dataframe if meta is provided
    if meta is not None:
        up_datasets = meta.loc[meta[meta_pertub_direction_col_name] == 'Up', meta_data_col_name].unique()
        down_datasets = meta.loc[meta[meta_pertub_direction_col_name] == 'Down', meta_data_col_name].unique()
    #run analysis for each data column in batch
    data_columns = [col for col in experiment.columns if col.startswith('data:')]
    res_batch = pd.DataFrame()
    for col in data_columns:
        #select data for column
        data = experiment[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE', col]]
        #determine enrichment type
        if meta is not None:
            if col in up_datasets:
                enrichment_type = 'enriched'
            elif col in down_datasets:
                enrichment_type = 'depleted'
            else:
                enrichment_type = enrichment_type
        #run analysis
        klRes = run_KL_diff_phos_enrichment_analysis(experiment=data, kin_type=kin_type, threshold=threshold, agg=agg, reformat_peptides=reformat_peptides, quant_type=quant_type, enrichment_type=enrichment_type, rename_kinases=rename_kinases, kinaseNames=kinaseNames)
        res_batch = pd.concat([res_batch, klRes], ignore_index=True)
    
    return res_batch
    


def run_KL_binary_enrichment_analysis(experiment, kin_type,
                                    determine_binary_foreground=True, reformat_peptides=True,
                                    meta=None, agg='mean', 
                                    threshold=1, evidence_size=None, 
                                    meta_data_col_name='benchmarking_column', meta_pertub_direction_col_name='perturbation_direction', 
                                    rename_kinases=False, kinaseNames=None):
    '''Run KL binary enrichment analysis on dataset by providing a list of foreground sites and
    using KL's general background phosphoproteome. For single data column.
    PARAMETERS:
    experiment: 
        dataframe with phosphoproteomic data
            data must be named in format 'data:dataset_name'
                should contain a single data column
            must contain columns for KSTAR_ACCESSION and KSTAR_SITE
            (if determine_binary_foreground=True)
    kin_type: str
        type of kinases to analyze ('tyrosine' or 'ser_thr')
    meta: dataframe with metadata for benchmarking dataset, must contain columns for dataset name and perturbation direction
    determine_binary_foreground: bool
        whether to determine foreground using Sam's benchmarking functions 
        if False: will to be provided with list of foreground sites
        (default: True) 
    reformat_peptides: bool
        whether to reformat peptides to central format for KL enrichment analysis (default: True)
    agg : {'count', 'mean'}
        method to use when aggregating duplicate substrate-sites. 
        'count' combines multiple representations and adds if values are non-NaN
        'mean' uses the mean value of numerical data from multiple representations of the same peptide.
            NA values are droped from consideration.
    threshold: numeric
        should be same quant type as data columns
            e.g. if data columns are fold change, threshold should be fold change value
        default is 1
    evidence_size: int or None
        number of sites to use as evidence for each dataset. If None, will use all sites that meet criteria for foreground.
        default is None
    meta_data_col_name: str
        column name in meta dataframe that contains column names in benchmarking dataset (should match the dataset names in the experiment dataframe)   
        default is 'benchmarking_column'    
    meta_pertub_direction_col_name: str
        column name in meta dataframe that contains expected direction of kinase activity change ('Up' or 'Down')
        default is 'perturbation_direction'
    rename_kinases: bool
        whether to rename kinases to standard names using kinaseNames resource file
    kinaseNames: dataframe with kinase names for mapping KL kinase names -> std names, must contain columns 'KL_name' and 'std_name' (required if rename_kinases=True)
    RETURNS:
    klRes: dataframe with KL enrichment results for data
    '''

    #DETERMINE FOREGROUND SITES
    if determine_binary_foreground:
        fg = benchmarking_functions.create_binary_evidence_for_benchmarking(experiment=experiment, meta=meta, agg=agg, 
                                                                            threshold=threshold, evidence_size=evidence_size, meta_data_col_name=meta_data_col_name, meta_pertub_direction_col_name=meta_pertub_direction_col_name)
        #add seq site info
        seq_info = experiment[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE']].drop_duplicates()
        fg = fg.merge(seq_info, on=['KSTAR_ACCESSION', 'KSTAR_SITE'], how='left')
        #get sites
        data_col = [col for col in experiment.columns if col.startswith('data:')][0]
        fg = pd.DataFrame(fg[fg[data_col] == True]['KSTAR_PEPTIDE'].unique(), columns=['Sequence'])
    else:
        #fg is a df with KSTAR_PEPTIDE column and column data:dataset_name with boolean values indicating whether site is in foreground or not
        data_col = [col for col in experiment.columns if col.startswith('data:')][0]
        fg = pd.DataFrame(experiment[experiment[data_col] == True]['KSTAR_PEPTIDE'].unique(), columns=['Sequence'])

    #REFORMAT PEPTIDES
    if reformat_peptides:
        #convert to central format 
        fg['Sequence'] = fg['Sequence'].apply(lambda x: convert_pep_to_central_format(x))
        #capitalize peptides
        fg['Sequence'] = fg['Sequence'].str.upper()

    #CREATE KL ENRICHMENT DATA OBJECT
    #drop any duplicate peptides in foreground
    fg = fg.drop_duplicates().reset_index(drop=True)
    enrich_data = kl.EnrichmentData(foreground=fg)

    #RUN KL ENRICHMENT ANALYSIS
    #determine kl_thresh based on kin_type
    if kin_type == 'tyrosine':
        kl_thresh = 8
    elif kin_type == 'ser_thr':
        kl_thresh = 15
    klRes = enrich_data.kinase_enrichment(kin_type=kin_type, kl_method='percentile_rank', kl_thresh=kl_thresh, non_canonical=True)
    #pull enrichment results dataframe from klRes object
    klRes = klRes.enrichment_results
    #add column for data column
    klRes['data_column'] = data_col
    #move index to column and name 'Kinase'
    klRes = klRes.reset_index().rename(columns={'index': 'kinase'})
    #rename kinases in klRes to standard names using kinaseNames resource file if rename_kinases=True
    if rename_kinases:
        if kinaseNames is None:
            raise ValueError('kinaseNames dataframe must be provided if rename_kinases=True')
        #make name mapping dict
        name_map = dict(zip(kinaseNames['KLETcst_py'], kinaseNames['ReName']))
        #store original names
        klRes['kinase_KLname'] = klRes['kinase'].copy()
        #map kinase names in klRes to standard names using name_map dict
        klRes['kinase'] = klRes['kinase_KLname'].map(name_map)

    return klRes


def run_KL_binary_enrichment_analysis_batch(experiment, kin_type, 
                                            determine_binary_foreground=True, meta=None,
                                            reformat_peptides=True,
                                            rename_kinases=False, kinaseNames=None,
                                             agg='mean', threshold=1, evidence_size=None,
                                             meta_data_col_name='benchmarking_column', meta_pertub_direction_col_name='perturbation_direction'):
    '''Run KL binary enrichment analysis on benchmarking dataset by providing a df of foreground sites and
    using KL's general background phosphoproteome. For multiple data columns in batch. 
    Parameters:
    experiment: dataframe with phosphoproteomic data
        data columns must be named in format 'data:dataset_name'
        must contain columns for KSTAR_ACCESSION and KSTAR_SITE
    kin_type: str
        type of kinase to analyze ('tyrosine' or 'ser_thr')
    determine_binary_foreground: bool
        whether to determine foreground using Sam's benchmarking functions 
        if False: will to be provided with list of foreground sites
        (default: True)
    reformat_peptides: bool
        whether to reformat peptides to central format for KL enrichment analysis (default: True)
    rename_kinases: bool
        whether to rename kinases in KL results to standard kinase names using kinaseNames resource file (default: False)
    kinaseNames: dataframe with kinase names for mapping KL kinase names -> std names, must contain columns 'KL_name' and 'std_name' (required if rename_kinases=True)
    meta: dataframe with metadata for benchmarking dataset, must contain columns for dataset name and perturbation direction
    agg : {'count', 'mean'}
        method to use when aggregating duplicate substrate-sites. 
        'count' combines multiple representations and adds if values are non-NaN
        'mean' uses the mean value of numerical data from multiple representations of the same peptide.
            NA values are droped from consideration.
        default is 'mean'
    threshold: numeric
        should be same quant type as data columns
            e.g. if data columns are fold change, threshold should be fold change value
        default is 1
    evidence_size: int or None
        number of sites to use as evidence for each dataset. If None, will use all sites that meet criteria for foreground.
        default is None
    meta_data_col_name: str
        column name in meta dataframe that contains column names in benchmarking dataset (should match the dataset names in the experiment dataframe)   
        default is 'benchmarking_column'
    meta_pertub_direction_col_name: str
        column name in meta dataframe that contains expected direction of kinase activity change ('Up' or 'Down')
        default is 'perturbation_direction'
    Returns:
    klRes: dataframe with KL enrichment results for each dataset in benchmarking data
    fg: dataframe with foreground sites
    '''
    
    #DETERMINE FOREGROUND SITES
    if determine_binary_foreground:
        fg = benchmarking_functions.create_binary_evidence_for_benchmarking(experiment=experiment, meta=meta, agg=agg, threshold=threshold, evidence_size=evidence_size, meta_data_col_name=meta_data_col_name, meta_pertub_direction_col_name=meta_pertub_direction_col_name)
        #add seq site info
        seq_info = experiment[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE']].drop_duplicates()
        fg = fg.merge(seq_info, on=['KSTAR_ACCESSION', 'KSTAR_SITE'], how='left')
    else:
        #fg is assumed to already be a binary df with unique rows for each site and columns for KSTAR_ACCESSION, KSTAR_SITE, and KSTAR_PEPTIDE
        fg = pd.DataFrame(experiment)

    #REFORMAT PEPTIDES
    if reformat_peptides:
        #convert to central format 
        fg['KSTAR_PEPTIDE'] = fg['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x))
        #capitalize peptides
        fg['KSTAR_PEPTIDE'] = fg['KSTAR_PEPTIDE'].str.upper()

    #RUN KL ENRICHMENT ANALYSIS FOR EACH DATA COLUMN
    klres = pd.DataFrame()
    data_columns = [col for col in experiment.columns if col.startswith('data:')]
    for data_col in data_columns:
        #select data for column
        data = fg[fg[data_col]][['KSTAR_PEPTIDE', data_col]]
        colRes = run_KL_binary_enrichment_analysis(experiment=data, kin_type=kin_type, 
                                                   determine_binary_foreground=False, reformat_peptides=False, meta=None, agg=agg, threshold=threshold, evidence_size=evidence_size, meta_data_col_name=meta_data_col_name, meta_pertub_direction_col_name=meta_pertub_direction_col_name)
        klres = pd.concat([klres, colRes], ignore_index=True)

    #RENAME KINASES IN KL RESULTS TO STANDARD NAMES
    if rename_kinases:
        if kinaseNames is None:
            raise ValueError('kinaseNames dataframe must be provided if rename_kinases=True')
        #make name mapping dict
        name_map = dict(zip(kinaseNames['KLETcst_py'], kinaseNames['ReName']))
        #store original names
        klres['kinase_KLname'] = klres['kinase'].copy()
        #map kinase names in klres to standard names using name_map dict
        klres['kinase'] = klres['kinase_KLname'].map(name_map)

    return klres

def run_KL_MEA_enrichment_analysis(experiment, kin_type, 
                                   agg='mean', reformat_peptides=True, rename_kinases=False, kinaseNames=None):
    '''Run KL motif enrichment analysis by providing df of phosphorylation data with peptide sequences and quantitative values. For single column of quantitative data
    PARAMETERS:
    experiment: dataframe with phosphoproteomic data
        must contain columns for KSTAR_ACCESSION, KSTAR_SITE, KSTAR_PEPTIDE, and at least one data column with quantitative values (e.g. fold change)
    kin_type: str
        type of kinase to analyze ('tyrosine' or 'ser_thr')
    agg: str
        aggregation method for duplicate sites ('mean', 'sum', False) (default: 'mean')
    reformat_peptides: bool
        whether to reformat peptides to central format for KL enrichment analysis (default: True)
    rename_kinases: bool
        whether to rename kinases in KL results to standard kinase names using kinaseNames resource file (default: False)   
    kinaseNames: dataframe with kinase names for mapping KL kinase names -> std names, must contain columns 'KL_name' and 'std_name' (required if rename_kinases=True)
    RETURNS:
    klRes: dataframe with KL enrichment results for data
    '''
    #REFORMAT EXPERIMENT
    #aggregate
    #if there are duplicate sites, will aggregate by mean value of data columns
    data_col = [col for col in experiment.columns if col.startswith('data:')][0]
    if agg:
        #aggregate data by peptide seq since KL MEA won't consider protein accession
        experiment = experiment.groupby([config.KSTAR_PEPTIDE])[data_col].agg(agg).reset_index()

    #reformat peptides to central format if reformat_peptides=True
    if reformat_peptides:
        experiment['KSTAR_PEPTIDE'] = experiment['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x))
        #capitalize peptides
        experiment['KSTAR_PEPTIDE'] = experiment['KSTAR_PEPTIDE'].str.upper()   
    
    #CREATE KL RANKED PHOSPHORYLATION DATA OBJECT
    #drop nans in data_col
    experiment = experiment.dropna(subset=[data_col])
    ranked_phos = kl.RankedPhosData(dp_data=experiment, rank_col=data_col, seq_col='KSTAR_PEPTIDE')

    #RUN KL ENRICHMENT ANALYSIS
    #determine kl_thresh based on kin_type
    if kin_type == 'tyrosine':
        kl_thresh = 8
    elif kin_type == 'ser_thr':
        kl_thresh = 15
    klRes = ranked_phos.mea(kin_type=kin_type, kl_method='percentile_rank', kl_thresh=kl_thresh, non_canonical=True)

    #pull enrichment results dataframe from klRes object
    klRes = klRes.enrichment_results
    #add column for data column
    klRes['data_column'] = data_col
    #move index to column and name 'Kinase'
    klRes = klRes.reset_index().rename(columns={'Kinase': 'kinase'})


    #RENAME KINASES IN KL RESULTS TO STANDARD NAMES
    if rename_kinases:
        if kinaseNames is None:
            raise ValueError('kinaseNames dataframe must be provided if rename_kinases=True')
        #make name mapping dict
        name_map = dict(zip(kinaseNames['KLETcst_py'], kinaseNames['ReName']))
        #store original names
        klRes['kinase_KLname'] = klRes['kinase'].copy()
        #map kinase names in klRes to standard names using name_map dict
        klRes['kinase'] = klRes['kinase_KLname'].map(name_map)

    return klRes


def run_KL_MEA_enrichment_analysis_batch(experiment, kin_type, 
                                         agg='mean', reformat_peptides=True, rename_kinases=False, kinaseNames=None):
    '''Run KL motif enrichment analysis by providing df of phosphorylation data with peptide sequences and quantitative values. For multiple columns of quantitative data in batch.
    PARAMETERS:
    experiment: dataframe with phosphoproteomic data
        must contain columns for KSTAR_ACCESSION, KSTAR_SITE, KSTAR_PEPTIDE, and at least one data column with quantitative values (e.g. fold change)
    kin_type: str
        type of kinase to analyze ('tyrosine' or 'ser_thr')
    agg: str
        aggregation method for duplicate sites ('mean', 'sum') (default: 'mean')
    reformat_peptides: bool
        whether to reformat peptides to central format for KL enrichment analysis (default: True)
    rename_kinases: bool
        whether to rename kinases in KL results to standard kinase names using kinaseNames resource file (default: False)   
    kinaseNames: dataframe with kinase names for mapping KL kinase names -> std names, must contain columns 'KL_name' and 'std_name' (required if rename_kinases=True)
    RETURNS:
    klRes: dataframe with KL enrichment results for each data column in experiment
    '''
    
    #REFORMAT EXPERIMENT
    #aggregate - by peptide seq since KL MEA won't consider protein accession
    data_columns = [col for col in experiment.columns if col.startswith('data:')]
    experiment = experiment.groupby([config.KSTAR_PEPTIDE])[data_columns].agg(agg).reset_index()
    #add KSTAR_PEPTIDE column back to experiment
    # seq_info = experiment[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE']].drop_duplicates(subset=['KSTAR_ACCESSION', 'KSTAR_SITE'])
    # experiment = experiment_agg.merge(seq_info, on=['KSTAR_ACCESSION', 'KSTAR_SITE'], how='left')

    #reformat peptides to central format if reformat_peptides=True
    if reformat_peptides:
        experiment['KSTAR_PEPTIDE'] = experiment['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x))
        #capitalize peptides
        experiment['KSTAR_PEPTIDE'] = experiment['KSTAR_PEPTIDE'].str.upper()   

    #RUN KL ENRICHMENT ANALYSIS FOR EACH DATA COLUMN
    klres = pd.DataFrame()
    for col in data_columns:
        #select data for column
        data = experiment[['KSTAR_PEPTIDE', col]] 
        colRes = run_KL_MEA_enrichment_analysis(experiment=data, kin_type=kin_type, agg=False, reformat_peptides=False, rename_kinases=False, kinaseNames=None)
        klres = pd.concat([klres, colRes], ignore_index=True)

    #REFORMAT RESULTS
    #rename kinases in klres to standard names using kinaseNames resource file if rename_kinases=True
    if rename_kinases:
        if kinaseNames is None:
            raise ValueError('kinaseNames dataframe must be provided if rename_kinases=True')
        #make name mapping dict
        name_map = dict(zip(kinaseNames['KLETcst_py'], kinaseNames['ReName']))
        #store original names
        klres['kinase_KLname'] = klres['kinase'].copy()
        #map kinase names in klres to standard names using name_map dict
        klres['kinase'] = klres['kinase_KLname'].map(name_map).fillna(klres['kinase_KLname'])

    return klres

def reformat_KL_fisher_enrichment_results(klRes):
    '''Reformat KL fisher enrichment results to KSTAR activities/sig results. 
    PARAMETERS:
    klRes: dataframe with KL differential phosphorylation enrichment results, from benchmarking functions using KL python functions
        should include column enrichment_type with values 'Up' or 'Down' to indicate direction of enrichment for each dataset
    RETURNS:
    kl_dp_activities: dataframe with log2 frequency factor from contigency table
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    kl_dp_sig: dataframe with FDR-corrected p-values from Fisher's exact test on contingency table
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    '''
    #REFORMAT KL RESULTS
    #rename kinase column to KSTAR_KINASE
    klRes = klRes.copy()
    klRes = klRes.rename(columns={'kinase': 'KSTAR_KINASE'})

    #REFORMAT RESULTS TO MATCH KSTAR ACTIVITIES/SIG FORMAT
    #store activities and sig in separate dataframes
    kl_dp_activities = pd.DataFrame()
    kl_dp_sig = pd.DataFrame()
    #reformat data for up and down regulated datasets separately, then combine at end
    for direction in ['Up', 'Down']:
        klRes_dir = klRes[klRes['enrichment_type'] == direction]
        #determine activities and sig columns to pivot based on enrichment_type col
        if direction == 'Up':
            activities_col = 'log2_freq_factor_upreg'
            sig_col = 'fisher_adj_pval_upreg'
        elif direction == 'Down':
            activities_col = 'log2_freq_factor_downreg'
            sig_col = 'fisher_adj_pval_downreg'
        #pivot to have kinases as index and data columns as columns, values are log2 frequency factor for activities and FDR-corrected p-value for sig
        activities = klRes_dir.pivot(index='KSTAR_KINASE', columns='data_column', values=activities_col)
        sig = klRes_dir.pivot(index='KSTAR_KINASE', columns='data_column', values=sig_col)
        #combine up and down regulated results
        kl_dp_activities = pd.concat([kl_dp_activities, activities], axis=1)
        kl_dp_sig = pd.concat([kl_dp_sig, sig], axis=1)

    return kl_dp_activities, kl_dp_sig

def reformat_KL_binary_enrichment_results(klRes):
    '''Reformat KL binary enrichment results to KSTAR activities/sig results.
    PARAMETERS:
    klRes: dataframe with KL binary enrichment results, from benchmarking functions using KL python functions
    RETURNS:
    kl_binary_activities: dataframe with log2 frequency factor from contigency table
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    kl_binary_sig: dataframe with FDR-corrected p-values from Fisher's exact test on contingency table
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    '''
    #REFORMAT KL RESULTS
    #rename kinase column to KSTAR_KINASE
    klRes = klRes.copy()
    klRes = klRes.rename(columns={'kinase': 'KSTAR_KINASE'})

    #PIVOT TO HAVE KINASES AS INDEX AND DATA COLUMNS AS COLUMNS, VALUES ARE LOG2 FREQUENCY FACTOR FOR ACTIVITIES AND FDR-CORRECTED P-VALUES FOR SIG
    kl_binary_activities = klRes.pivot(index='KSTAR_KINASE', columns='data_column', values='log2_freq_factor')
    kl_binary_sig = klRes.pivot(index='KSTAR_KINASE', columns='data_column', values='fisher_adj_pval')

    return kl_binary_activities, kl_binary_sig

def reformat_KL_MEA_enrichment_results(klRes):
    '''Reformat KL MEA enrichment results to KSTAR activities/sig results.
    PARAMETERS:
    klRes: dataframe with KL MEA enrichment results, from benchmarking functions using KL python functions
    RETURNS:
    kl_MEA_activities: dataframe with log2 frequency factor from contigency table
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    kl_MEA_sig: dataframe with FDR-corrected p-values from Fisher's exact test on contingency table
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    '''
    #REFORMAT KL RESULTS
    #rename kinase column to KSTAR_KINASE
    klRes = klRes.copy()
    klRes = klRes.rename(columns={'kinase': 'KSTAR_KINASE'})

    #PIVOT TO HAVE KINASES AS INDEX AND DATA COLUMNS AS COLUMNS, VALUES ARE LOG2 FREQUENCY FACTOR FOR ACTIVITIES AND FDR-CORRECTED P-VALUES FOR SIG
    kl_MEA_activities = klRes.pivot(index='KSTAR_KINASE', columns='data_column', values='NES')
    kl_MEA_sig = klRes.pivot(index='KSTAR_KINASE', columns='data_column', values='FDR')

    return kl_MEA_activities, kl_MEA_sig

def reformat_KL_res(klRes):
    '''Reformat KL enrichment results to KSTAR activities/sig results. Determines type of KL enrichment analysis based on columns in klRes and reformats accordingly.
    PARAMETERS:
    klRes: dataframe with KL enrichment results, from benchmarking functions using KL python functions
    RETURNS:
    kl_activities: dataframe with log2 frequency factor from contigency table or NES for MEA
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    kl_sig: dataframe with FDR-corrected p-values from Fisher's exact test on contingency table or MEA FDR
            index: KSTAR_KINASE
            cols: data_col_1, data_col_2
    '''
    if 'log2_freq_factor_upreg' in klRes.columns and 'log2_freq_factor_downreg' in klRes.columns:
        return reformat_KL_fisher_enrichment_results(klRes)
    elif 'log2_freq_factor' in klRes.columns:
        return reformat_KL_binary_enrichment_results(klRes)
    elif 'NES' in klRes.columns and 'FDR' in klRes.columns:
        return reformat_KL_MEA_enrichment_results(klRes)
    else:
        raise ValueError('klRes does not contain expected columns for KL enrichment results. Please check columns and ensure klRes is in expected format.')

#%%LOAD DATA AND RESOURCES

#load benchmarking data from KSTAR paper
dataDir = './Data/Perturbation'
#Y
path = os.path.join(dataDir, 'benchmarking_dataset_Y.tsv')
yData = pd.read_csv(path, sep='\t')
#ST
path = os.path.join(dataDir, 'benchmarking_dataset_ST.tsv')
stData = pd.read_csv(path, sep='\t')

#load metadata for benchmarking data
#Y
path = os.path.join(dataDir, 'meta_y.csv')
metaY = pd.read_csv(path)
#ST
path = os.path.join(dataDir, 'meta_st.csv')
metaST = pd.read_csv(path)

#load kinase names for mapping KL kinase names -> std names
resourceDir = './Data/Resources'
path = os.path.join(resourceDir, 'KinaseNames_KLPSSMs_NK_KLETmit_KLETcstMEA.csv')
kinaseNames = pd.read_csv(path)

#%%MAKE RESULTS DIRECTORY
resDir = 'Analysis_Results'
if not os.path.exists(resDir):
    os.makedirs(resDir)

#%%KL DIFFERENTIAL PHOSPHORYLATION ENRICHMENT ANALYSIS #################################
#%%make dir to save KL diff phos enrichment results
resSubDir = 'KL_DiffPhos_Enrichment_Results'
diffPhosResDir = os.path.join(resDir, resSubDir)
if not os.path.exists(diffPhosResDir):
    os.makedirs(diffPhosResDir)

#%%run KL diff phos enrichment analysis for Y data
klDiffPhosResY = run_KL_diff_phos_enrichment_analysis_batch(experiment=yData, kin_type='tyrosine', threshold=1, agg='mean', reformat_peptides=True, quant_type='fold_change', enrichment_type='both', meta=metaY, meta_data_col_name='benchmarking_column', meta_pertub_direction_col_name='perturbation_direction', rename_kinases=True, kinaseNames=kinaseNames)

klDiffPhosResY.to_csv(os.path.join(diffPhosResDir, 'KL_diff_phos_enrichment_results_Y.csv'), index=False)

print ('Analysis complete! Results saved to:', diffPhosResDir)

#%%run KL diff phos enrichment analysis for ST data
klDiffPhosResST = run_KL_diff_phos_enrichment_analysis_batch(experiment=stData, kin_type='ser_thr', threshold=1, agg='mean', reformat_peptides=True, quant_type='fold_change', enrichment_type='both', meta=metaST, meta_data_col_name='benchmarking_column', meta_pertub_direction_col_name='perturbation_direction', rename_kinases=True, kinaseNames=kinaseNames)

klDiffPhosResST.to_csv(os.path.join(diffPhosResDir, 'KL_diff_phos_enrichment_results_ST.csv'), index=False)

print('Analysis complete! Results saved to:', diffPhosResDir)

#%%KL BINARY ENRICHMENT ANALYSIS ###############################################
#%%#make dir to save KLenrichment results
resSubDir = 'KL_Enrichment_Results'
enrResDir = os.path.join(resDir, resSubDir)
if not os.path.exists(enrResDir):
    os.makedirs(enrResDir)

#%%run KL enrichment analysis for Y data
klResY = run_KL_binary_enrichment_analysis_batch(experiment=yData, kin_type='tyrosine', determine_binary_foreground=True, meta=metaY, reformat_peptides=True, rename_kinases=True, kinaseNames=kinaseNames)
klResY.to_csv(os.path.join(enrResDir, 'KL_enrichment_results_Y.csv'), index=False)

print ('Analysis complete! Results saved to:', enrResDir)

#%%run KL enrichment analysis for ST data
klResST = run_KL_binary_enrichment_analysis_batch(experiment=stData, kin_type='ser_thr', determine_binary_foreground=True, meta=metaST, reformat_peptides=True, rename_kinases=True, kinaseNames=kinaseNames)
klResST.to_csv(os.path.join(enrResDir, 'KL_enrichment_results_ST.csv'), index=False)

print('Analysis complete! Results saved to:', enrResDir)

#%%MEA: MOTIF ENRICHMENT ANALYSIS ###############################################
#make dir to save MEA results
resSubDir = 'KL_MEA_Enrichment_Results'
enrResDir = os.path.join(resDir, resSubDir)
if not os.path.exists(enrResDir):
    os.makedirs(enrResDir)
#%%run KL motif enrichment analysis for Y data
MEAresY = run_KL_MEA_enrichment_analysis_batch(experiment=yData, kin_type='tyrosine', agg='mean', reformat_peptides=True, rename_kinases=True, kinaseNames=kinaseNames)
MEAresY.to_csv(os.path.join(enrResDir, 'KL_MEA_enrichment_results_Y.csv'), index=False)

print ('Analysis complete! Results saved to:', enrResDir)

#%%run KL motif enrichment analysis for ST data
MEAresST = run_KL_MEA_enrichment_analysis_batch(experiment=stData, kin_type='ser_thr', agg='mean', reformat_peptides=True, rename_kinases=True, kinaseNames=kinaseNames)
MEAresST.to_csv(os.path.join(enrResDir, 'KL_MEA_enrichment_results_ST.csv'), index=False)

print('Analysis complete! Results saved to:', enrResDir)

# %%MAKE SINGLE COLUMN VERSION OF BENCHMARKING DATA TO DOUBLE CHECK VS KL WEB APP RESULTS
#make dir to save reformatted benchmarking data
resSubDir = 'Reformatted_Benchmarking_Data'
reformatDataDir = os.path.join(resDir, resSubDir)
if not os.path.exists(reformatDataDir):
    os.makedirs(reformatDataDir)

#FOR DIFFERENTIAL PHOSPHORYLATION
#reformat Y data
#get information columns and single test column
info_cols = [col for col in yData.columns if col.startswith('KSTAR_')]
test_col = [col for col in yData.columns if col.startswith('data:')][0]
yData_reformat = yData[info_cols + [test_col]]
#drop nas in test_col
yData_reformat = yData_reformat.dropna(subset=[test_col])
#aggregate by mean if there are duplicate sites
yData_reformat = yData_reformat.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[test_col].agg('mean').reset_index()
#add small value to test_col to avoid issues with log2 transformation of 0 values
yData_reformat[test_col] = yData_reformat[test_col].replace(0, 1e-10)
yData_reformat[test_col] = np.log2(yData_reformat[test_col])
#add KSTAR_PEPTIDE
peptides = yData[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE']].drop_duplicates()
yData_reformat = yData_reformat.merge(peptides, on=['KSTAR_ACCESSION', 'KSTAR_SITE'])
#reformat peptides to central format and capitalize
yData_reformat['KSTAR_PEPTIDE'] = yData_reformat['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x)).str.upper()
#drop rows with missing values in test_col
yData_reformat = yData_reformat.dropna(subset=[test_col]).reset_index(drop=True)
#save reformatted data
yData_reformat.to_csv(os.path.join(reformatDataDir, f'Y_data_reformatted_DP_{test_col}.tsv'), index=False, sep='\t')

#FOR BINARY ENRICHMENT
#reformat Y data
#get information columns and single test column
info_cols = [col for col in yData.columns if col.startswith('KSTAR_')]
test_col = [col for col in yData.columns if col.startswith('data:')][0]
yData_binary_reformat = yData[info_cols + [test_col]]
#drop nas in test_col
yData_binary_reformat = yData_binary_reformat.dropna(subset=[test_col])
#make binary column for test_col based on threshold of 1 (can adjust if needed)
yData_binary_reformat = benchmarking_functions.create_binary_evidence_for_benchmarking(experiment=yData_binary_reformat, meta=metaY, agg='mean', threshold=1, evidence_size=None, meta_data_col_name='benchmarking_column', meta_pertub_direction_col_name='perturbation_direction')
#add KSTAR_PEPTIDE
peptides = yData[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_PEPTIDE']].drop_duplicates()
yData_binary_reformat = yData_binary_reformat.merge(peptides, on=['KSTAR_ACCESSION', 'KSTAR_SITE'], how='left')
#reformat peptides to central format and capitalize
yData_binary_reformat['KSTAR_PEPTIDE'] = yData_binary_reformat['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x)).str.upper()
#keep only peptides where binary column is True, no duplicates
yData_binary_reformat = yData_binary_reformat[yData_binary_reformat[test_col] == True][['KSTAR_PEPTIDE']].drop_duplicates().reset_index(drop=True)
#save reformatted data 
yData_binary_reformat.to_csv(os.path.join(reformatDataDir, f'Y_data_reformatted_binary_{test_col}.tsv'), index=False, sep='\t')

#OCHOA Phosphoproteome - bg
#load Ochoa phosphoproteome data
path = '../Naegle_Lab/Networks/AnalysisResults/tyrosine_phosphoproteome_scored_Ochoa_reformatted.csv'
ochoa_bg = pd.read_csv(path)
#keep only columns for sequence and save as tsv
ochoa_bg = ochoa_bg[['pep_KL_SS_PSOPL']].drop_duplicates().dropna().reset_index(drop=True).rename(columns={'pep_KL_SS_PSOPL': 'Sequence'})
#reformat peptides to central format and capitalize
ochoa_bg['Sequence'] = ochoa_bg['Sequence'].apply(lambda x: convert_pep_to_central_format(x)).str.upper()
ochoa_bg.to_csv(os.path.join(reformatDataDir, 'Ochoa_phosphoproteome_background.tsv'), index=False, sep='\t')
ochoa_bg

#FOR MEA
#get information columns and single test column
info_cols = [col for col in yData.columns if col.startswith('KSTAR_')]
test_col = [col for col in yData.columns if col.startswith('data:')][0]
yData_mea_reformat = yData[info_cols + [test_col]]
#drop nas in test_col
yData_mea_reformat = yData_mea_reformat.dropna(subset=[test_col])
#aggregate by mean if there are duplicate sites
yData_mea_reformat = yData_mea_reformat.groupby([config.KSTAR_PEPTIDE])[test_col].agg('mean').reset_index()
# #convert test_col to log2 fold change
# yData_mea_reformat[test_col] = yData_mea_reformat[test_col].replace(0, 1e-10)
# yData_mea_reformat[test_col] = np.log2(yData_mea_reformat[test_col])
#reformat peptides to central format and capitalize
yData_mea_reformat['KSTAR_PEPTIDE'] = yData_mea_reformat['KSTAR_PEPTIDE'].apply(lambda x: convert_pep_to_central_format(x)).str.upper()
#save reformatted data
yData_mea_reformat.to_csv(os.path.join(reformatDataDir, f'Y_data_reformatted_MEA_{test_col}.tsv'), index=False, sep='\t')

#%%CALCULATE ACCURACY (PHIT) 
#%%#address meta data

#expand meta so that each row corresponds to a single kinase condition
y_meta = metaY.copy()
y_meta['Kinases'] = y_meta['kinases'].apply(lambda x: x.split(','))
y_meta = y_meta.explode('Kinases')
y_meta.rename({'benchmarking_column': 'Dataset'}, axis = 1, inplace = True)
y_meta = y_meta.drop_duplicates()
y_meta.reset_index(inplace = True)

#expand meta so that each row corresponds to a single kinase condition
st_meta = metaST.copy()
st_meta['Kinases'] = st_meta['kinases'].apply(lambda x: str(x).split(','))
st_meta = st_meta.explode('Kinases')
st_meta.rename({'benchmarking_column': 'Dataset'}, axis = 1, inplace = True)
st_meta = st_meta.drop_duplicates()
st_meta.reset_index(inplace = True)
meta = {'Y':y_meta,'ST':st_meta}

#number of tests
total_tests = meta['Y'].shape[0]+meta['ST'].shape[0]

#%%Organize results
#load activity scores/ranks and combine
#######################################
#### KL diff phospho results ####
path = os.path.join(diffPhosResDir, 'KL_diff_phos_enrichment_results_Y.csv')
KL_DP_y = pd.read_csv(path)
#reformat to KSTAR activities/sig format
KL_DP_y, KL_DP_y_sig = reformat_KL_fisher_enrichment_results(KL_DP_y)

path = os.path.join(diffPhosResDir, 'KL_diff_phos_enrichment_results_ST.csv')
KL_DP_st = pd.read_csv(path)
#reformat to KSTAR activities/sig format
KL_DP_st, KL_DP_st_sig = reformat_KL_fisher_enrichment_results(KL_DP_st)

# activities
KL_DP_activities = {'Y':KL_DP_y,'ST':KL_DP_st}
#sig 
KL_DP_sig = {'Y':KL_DP_y_sig, 'ST': KL_DP_st_sig}

#### KL binary enrichment results ####
path = os.path.join(enrResDir, 'KL_enrichment_results_Y.csv')
KL_BE_y = pd.read_csv(path)
#reformat to KSTAR activities/sig format
KL_BE_y, KL_BE_y_sig = reformat_KL_binary_enrichment_results(KL_BE_y)

path = os.path.join(enrResDir, 'KL_enrichment_results_ST.csv')
KL_BE_st = pd.read_csv(path)
KL_BE_st, KL_BE_st_sig = reformat_KL_binary_enrichment_results(KL_BE_st)

#activities
KL_BE_activities = {'Y':KL_BE_y, 'ST': KL_BE_st}
#sig
KL_BE_sig = {'Y':KL_BE_y_sig, 'ST': KL_BE_st_sig}

#### KL MEA results ####
path = os.path.join(enrResDir, 'KL_MEA_enrichment_results_Y.csv')
KL_MEA_y = pd.read_csv(path)
#reformat to KSTAR activities/sig format
KL_MEA_y, KL_MEA_y_sig = reformat_KL_MEA_enrichment_results(KL_MEA_y)

path = os.path.join(enrResDir, 'KL_MEA_enrichment_results_ST.csv')
KL_MEA_st = pd.read_csv(path)
KL_MEA_st, KL_MEA_st_sig = reformat_KL_MEA_enrichment_results(KL_MEA_st)

#activities
KL_MEA_activities = {'Y':KL_MEA_y, 'ST': KL_MEA_st}
#sig
KL_MEA_sig = {'Y':KL_MEA_y_sig, 'ST': KL_MEA_st_sig}

#### combine into a single dictionary (for if you are comparing multiple results)
activities = {'DiffPhos':KL_DP_activities, 'Binary':KL_BE_activities, 'MEA':KL_MEA_activities}
sig = {'DiffPhos':KL_DP_sig, 'Binary':KL_BE_sig, 'MEA':KL_MEA_sig}

# %%##### Save data to combine with other benchmarking results 
#make dir to save reformatted benchmarking results
phit_dictionaries = 'Phit_dictionaries'
phitDir = os.path.join(resDir, phit_dictionaries)
if not os.path.exists(phitDir):
    os.makedirs(phitDir)

#save activities and sig dictionaries as pickle files
with open(os.path.join(phitDir, 'activities_dict.p'), 'wb') as f:
    pickle.dump(activities, f)
with open(os.path.join(phitDir, 'sig_dict.p'), 'wb') as f:
    pickle.dump(sig, f)
#%% Get hits matrix
#get hits matrix
hits_matrix_rank = benchmarking_functions.get_hit_matrix(meta, results = activities, result_type = 'rank', k = 10)
hits_matrix_sig, dataset_hits = benchmarking_functions.get_hit_matrix(meta, results = sig, result_type = 'sig', p = 0.05)

#%%Calculate accuracy (Phit)
#calculate Phit for each kinase
Phit_rank = benchmarking_functions.calculate_Phit(hits_matrix_rank)
Phit_sig = benchmarking_functions.calculate_Phit(hits_matrix_sig)



# %%Plot Phit results
benchmarking_functions.Phit_barplot(Phit_rank, Phit_sig)

# %%
Phit_rank
# %%
benchmarking_functions.kinase_specific_heatmap(hits_matrix_rank, hits_matrix_sig, mod = 'Y', figsize = (6,2))
# %%
