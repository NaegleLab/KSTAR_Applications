import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import matplotlib.pyplot as plt

class KSEA:
    def __init__(self, experiment, network, data_cols = None, kinase_col = 'Kinase', log_transform = True, threshold = None, threshold_on = None):
        #save network and experiment
        self.network = network.copy()
        self.experiment = experiment.copy()
            
        
        #Treshold the network if thresholding
        if threshold is not None:
            self.network = self.network[self.network[threshold_on] >= threshold]
            self.network.drop([threshold_on], axis = 1, inplace = True)
            print('Network thresholded')

       
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
        
        #If data_cols is none, look for columns starting with 'data:'
        if data_cols is None:
            self.data_cols = [col for col in self.experiment.columns if 'data:' in col]
        else:
            self.data_cols = data_cols
            
                #if data hasn't been log transformed, do so now
        if log_transform:
            self.experiment[self.data_cols] = np.log2(self.experiment[self.data_cols])
            #if an values returned -inf (0s in the dataframe), change to NaN and notify user
            if -np.inf in self.experiment[self.data_cols].values:
                self.experiment.replace([-np.inf, np.inf], np.nan, inplace = True)
                print('Log transform produced -inf values (0s in the data columns). Replaced with NaN.')
       

        self.all_results = {}
        self.zscore = pd.DataFrame(None, columns = self.data_cols, index = self.kinases)
        self.pvals = pd.DataFrame(None, columns = self.data_cols, index = self.kinases)
        self.FDR = pd.DataFrame(None, columns = self.data_cols, index = self.kinases)
        self.m = pd.DataFrame(None, columns = self.data_cols, index = self.kinases)
        
    def runKSEA_singleExperiment(self, data_col):
        """
        Perform KSEA analysis for each kinase with substrates identified in experiment, with quantitative information
        found in  data_col
        
        Parameters
        ----------
        data_col: string
            which data column to perform KSEA analysis on
            
        Returns
        -------
        all_results: pandas dataframe
            melted dataframe including all KSEA results, including zscore, number of identified substrates (m), and false discovery rate
        zscore: pandas dataframe
            dataframe which contains zscores from all data columns, updated with new zscores calculated for data_col
        pvals: pandas dataframe
            dataframe which contains p-values from all data columns, updated with new p-values calculated for data_col
        FDR: pandas dataframe
            dataframe which contains false discover rate from all data columns, updated with new false discovery rate calculated for data_col
        m: pandas dataframe
            dataframe which contains number of identified substrates for each kianse from all data columns, updated with new m calculated for data_col
        """
        #save experiment in a temporary df, dropna columns
        tmp_exp = self.experiment.dropna(subset = [data_col], axis = 0)
        #rename data column to log2FC
        tmp_exp.rename({data_col:'log2FC'}, inplace = True, axis = 1)
        #if experiment has the same kinase column as network, remove it (don't need it)
        if self.kinase_col in tmp_exp.columns:
            tmp_exp.drop(self.kinase_col, axis = 1, inplace = True)
        #merge network and experiment, keep important columns, and average log2FC for the same sites
        merged_network = self.network.merge(tmp_exp, on = ['KSTAR_ACCESSION', 'KSTAR_SITE'])
        merged_network = merged_network[[self.kinase_col,'KSTAR_ACCESSION', 'KSTAR_SITE', 'log2FC']]
        network_agg = merged_network.groupby(by = [self.kinase_col,'KSTAR_ACCESSION','KSTAR_SITE']).aggregate('mean').reset_index()
        
        #perform KSEA analysis, as implemented in Casado et al.
        meanFC = network_agg.groupby(by = self.kinase_col)['log2FC'].mean().reset_index()
        meanFC['Enrichment'] = meanFC['log2FC']/abs(np.mean(tmp_exp['log2FC']))
        meanFC['m'] = network_agg.groupby(by = self.kinase_col)['log2FC'].count().values
        meanFC['zscore'] = (meanFC['log2FC'] - np.mean(tmp_exp['log2FC']))*(meanFC['m'])**0.5/np.std(tmp_exp['log2FC'])
        meanFC['pval'] = stats.norm.cdf(-abs(meanFC['zscore'].values))
        meanFC['FDR'] = fdrcorrection(meanFC['pval'])[1]
        
        
        self.all_results[data_col] = meanFC
        self.zscore.loc[meanFC[self.kinase_col].values, data_col] = meanFC['zscore'].values
        self.pvals.loc[meanFC[self.kinase_col].values, data_col] = meanFC['pval'].values
        self.FDR.loc[meanFC[self.kinase_col].values, data_col] = meanFC['FDR'].values
        self.m.loc[meanFC[self.kinase_col].values, data_col] = meanFC['m'].values
        
    def runKSEA(self):
        """
        Perform KSEA analysis on all data columns in the experiment
        """
        for col in self.data_cols:
            self.runKSEA_singleExperiment(col)
        
        #drop kinases that don't have any predictions in any experiment
        self.zscore.dropna(how = 'all', inplace = True)
        self.pvals.dropna(how = 'all', inplace = True)
        self.FDR.dropna(how = 'all', inplace = True)
        self.m.dropna(how = 'all', inplace = True)
        
    def saveResults(self, fname = 'ksea', odir = ''):
        """
        Saves zscores, FDR, and m for all data columns in seperate .tsv files
        """
        self.zscore.to_csv(f'{odir}{fname}_zscore.tsv', sep = '\t')
        self.FDR.to_csv(f'{odir}{fname}_FDR.tsv', sep = '\t')
        self.m.to_csv(f'{odir}{fname}_m.tsv', sep = '\t')
        
        
        
def plotResults(KSEA, name = None, figsize = (10,10), ax = None, corrected = True):
    """
    Plots KSEA zscores for a single data column using a horizontal bar chart, sorted according to zscores
    and with significant bars colored with red
    """
    if name is None:
        print('Please indicate which result you would like to plot')
    else:
        cmap = sns.color_palette('colorblind')
        
        if ax is None:
            fig, ax = plt.subplots(figsize = figsize)
        
        plot_data = KSEA.zscore[name].sort_values(ascending = True)
        plot_data.dropna(inplace = True)
        kinases = plot_data.index

        #plot
        bar = ax.barh(kinases, plot_data, color = '#6B838F')
        
        #color Zscore with p <= 0.05
        if corrected:
            significance = (KSEA.FDR.loc[kinases,name] <= 0.05)
            sig_loc = np.where(significance)[0]
        else:
            significance = (KSEA.pvals.loc[kinases,name] <= 0.05)
            sig_loc = np.where(significance)[0]
            
        for loc in sig_loc:
            bar[loc].set(color = '#FF3300')

        xlabel = plt.xlabel('Z-score')