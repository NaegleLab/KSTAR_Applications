import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

class KARP:
    def __init__(self, experiment, network, data_cols = None, kinase_col = 'Kinase', log_transform = False, threshold_on = None):
        #save network and experiment
        self.network = network.copy()
        self.experiment = experiment.copy()
        
        
        #find overlap between network and experiment
        self.exp_substrates = list(experiment['KSTAR_ACCESSION']+'_'+experiment['KSTAR_SITE'])
        self.net_substrates = list(network['KSTAR_ACCESSION']+'_'+network['KSTAR_SITE'])
        self.site_overlap = list(set(self.exp_substrates).intersection(set(self.net_substrates)))
        
        self.kinase_col = kinase_col
        #find unique kinases
        self.kinases = network[kinase_col].unique()
        
        #extract substrate info for each kinase in the network
        self.ks_info = {}
        for kin in self.kinases:
            tmp_net = network[network[kinase_col] == kin]
            kinase_substrates = list(tmp_net['KSTAR_ACCESSION']+'_'+tmp_net['KSTAR_SITE'])
            self.ks_info[kin] = kinase_substrates
        
        #find total number of substrates for each kinase in PSP
        self.t = self.network.groupby(by = self.kinase_col).sum().reset_index()
        
        #If data_cols is none, look for columns starting with 'data:'
        if data_cols is None:
            self.data_cols = [col for col in self.experiment.columns if 'data:' in col]
        else:
            self.data_cols = data_cols
            
                #if data hasn't been log transformed, do so now
        #if log_transform:
        #    self.experiment[self.data_cols] = np.log2(self.experiment[self.data_cols])
        #    #if an values returned -inf (0s in the dataframe), change to NaN and notify user
        #    if -np.inf in self.experiment[self.data_cols].values:
        #        self.experiment.replace([-np.inf, np.inf], np.nan, inplace = True)
        #        print('Log transform produced -inf values (0s in the data columns). Replaced with NaN.')
       

        self.all_results = {}
        self.K = pd.DataFrame(None, columns = self.data_cols, index = self.kinases)
        self.rank = pd.DataFrame(None, columns = self.data_cols, index = self.kinases)
        
        
    def runKARP_singleExperiment(self, data_col):
        """
        Perform KARP analysis for each kinase with substrates identified in experiment, with quantitative information
        found in  data_col
        
        Parameters
        ----------
        data_col: string
            which data column to perform KSEA analysis on
            
        Returns
        -------
        all_results: pandas dataframe
            melted dataframe including all KARP results including K-score and rank
        K: pandas dataframe
            dataframe which contains K from all data columns, updated with new K-scores calculated for data_col
        rank: pandas dataframe
            dataframe which contains ramnk from all data columns, updated with new p-values calculated for data_col
      
        """
        #save experiment in a temporary df, dropna columns
        tmp_exp = self.experiment.dropna(subset = [data_col], axis = 0)
        #if experiment has the same kinase column as network, remove it (don't need it)
        if self.kinase_col in tmp_exp.columns:
            tmp_exp.drop(self.kinase_col, axis = 1, inplace = True)
        #merge network and experiment, keep important columns, and average log2FC for the same sites
        merged_network = self.network.merge(tmp_exp, on = ['KSTAR_ACCESSION', 'KSTAR_SITE'])
        merged_network = merged_network[[self.kinase_col,'KSTAR_ACCESSION', 'KSTAR_SITE', data_col]]
        network_agg = merged_network.groupby(by = [self.kinase_col,'KSTAR_ACCESSION','KSTAR_SITE']).aggregate('mean').reset_index()
        
        #perform KARP analysis
        
        #aggregate
        sum_substrates = network_agg.groupby(by = self.kinase_col)[data_col].sum()
        sum_phosphosites = tmp_exp[data_col].sum()
        m = network_agg.groupby(by = self.kinase_col).count()['KSTAR_SITE']
        t = self.network.groupby(by = self.kinase_col).count()['KSTAR_SITE']
        
        #calculate K-score
        self.K.loc[network_agg[self.kinase_col].values,data_col] = (sum_substrates/sum_phosphosites)*(m/t.loc[m.index])**0.5*10**6
        

        
        
        
        
        
        
    def runKARP(self):
        """
        Perform KARP analysis on all data columns in the experiment
        """
        for col in self.data_cols:
            self.runKARP_singleExperiment(col)
        

        self.K.dropna(how = 'all', inplace = True)
        
        
        
