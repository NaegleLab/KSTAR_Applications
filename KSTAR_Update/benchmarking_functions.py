import pandas as pd
import numpy as np

#plotting
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.tri import Triangulation
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patches as patches

#stats
import scipy.stats as stats
import math
import random

from kstar import config
from kstar.plot import DotPlot

def file_to_list(file):
    """
    Read a file and return a list for each entry

    Parameters
    ----------
    file : str
        File to read
    
    Returns
    -------
    lst : list
        List of entries in the file
    """
    with open(file, 'r') as f:
        lst = f.readlines()
    lst = [x.strip() for x in lst]
    return lst

def list_to_file(lst, file):
    """
    Write a list to a file

    Parameters
    ----------
    lst : list
        List to write to file
    file : str
        File to write to
    """
    with open(file, 'w') as f:
        for item in lst:
            f.write("%s\n" % item)
    
###########################################################################################################################################
######################## Accuracy Analysis #############################################################################################
###########################################################################################################################################

def create_binary_evidence_for_benchmarking(experiment, meta, agg = 'mean', threshold = 1.0,  evidence_size = None, 
                                            meta_data_col_name = 'Dataset', meta_pertub_direction_col_name = 'Direction'):
    """
    Returns a binary evidence data frame from the benchmarking dataset according to the parameters passed in for method for aggregating
    duplicates and considering whether a site is included as evidence or not. Adapted from kinase activity class function, but applies different thresholds based on if the data is inhibition or

    Parameters
    ----------
    threshold : float
        threshold value used to filter rows
    evidence_size: None or int
        the number of sites to use for prediction for each sample. If a value is provided, this will override the threshold, and will instead obtain the N sites with the greatest abundance within each sample.
    agg : {'count', 'mean'}
        method to use when aggregating duplicate substrate-sites. 
        'count' combines multiple representations and adds if values are non-NaN
        'mean' uses the mean value of numerical data from multiple representations of the same peptide.
            NA values are droped from consideration.
    meta_data_col_name: str
        column name in meta dataframe that contains column names in benchmarking dataset (should match the dataset names in the experiment dataframe)
    meta_pertub_direction_col_name: str
        column name in meta dataframe that contains expected direction of kinase activity change ('Up' or 'Down')

    Returns
    -------
    evidence_binary : pd.DataFrame
        Matches the evidence dataframe of the kinact object, but with 0 or 1 if a site is included or not.
        This is uniquified and rows that are never used are removed.
    
        
    """
    #check to make sure datasets have data
    data_columns = [col for col in experiment.columns if 'data:' in col]
    #check_data_columns(experiment, data_columns, agg, threshold, greater)

    #grab stimulation and inhibition datasets
    #split data into datasets with upregulated kinases and downregulated kinases
    up_datasets = meta.loc[meta[meta_pertub_direction_col_name] == 'Up', meta_data_col_name].unique()
    down_datasets = meta.loc[meta[meta_pertub_direction_col_name] == 'Down', meta_data_col_name].unique()

    #collapse sites into single row based on agg parameter
    experiment = experiment.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[data_columns].agg(agg).reset_index()
    
    #set the binary evidence for whether a site is included
    evidence_binary = experiment.copy()
    if evidence_size is None:
        for col in data_columns:
            if col in up_datasets:
                evidence_binary[col] = evidence_binary[col] >= threshold
            elif col in down_datasets:
                evidence_binary[col] = evidence_binary[col] <= threshold
    else:
        for col in data_columns:
            #check how many non_nan sites there (if less than N, set n to be equal to number of sites available
            num_sites_available = evidence_binary.dropna().shape[0]
            if num_sites_available >= evidence_size:
                n = evidence_size
            else:
                n = num_sites_available
                print(f"{col} has less than {evidence_size} sites available, so using {n} sites instead")
                
            if col in up_datasets:
                max_indices = np.argsort(-evidence_binary[col].values)[0:n]
                evidence_binary[col] = 0
                col_loc = np.where(evidence_binary.columns == col)[0][0]
                evidence_binary.iloc[max_indices, col_loc] = 1
            elif col in down_datasets:
                min_indices = np.argsort(evidence_binary[col].values)[0:n]
                evidence_binary[col] = 0
                col_loc = np.where(evidence_binary.columns == col)[0][0]
                evidence_binary.iloc[min_indices, col_loc] = 1

    #remove phosphorylation sites that were not selected in any experiment (useful for very large experiments where removing the need to copy data reduces time)
    evidence_binary.drop(evidence_binary[evidence_binary[data_columns].sum(axis=1) == 0].index, inplace = True) 

    #add back compendia/study bias information to binary evidence
    compendia = config.HUMAN_REF_COMPENDIA[['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS']]
    evidence_binary = evidence_binary.merge(compendia, on = ['KSTAR_ACCESSION', 'KSTAR_SITE'], how = 'left')
        
        
    return evidence_binary

def get_hit_matrix(meta, results = {}, result_type = 'rank', k = 10, p =0.05):
    """
    Function to generate matrix indicating whether or not KSTAR successfully predicted a kinase as differentially active in a given perturbation dataset from the global benchmarking data used in the main publication. 1 indicates a hit, 0 indicates a miss, and NaN indicates that the kinase did not have predictions for KSTAR.

    Parameters
    ----------
    results : nested dict
        Dictionary containing the results of different implementations of KSTAR, either the activity (-log10 of Mann Whitney) or FPR. If activity, result_type should be 'rank'. If FPR, result_type should be 'sig'. Should be in the format of results[algorithm][modification] = pd.DataFrame (in other words, when creating the dictionary, the first key should be the name of the algorithm, and the second key should be the modification type (Y or ST). The value should be a pandas dataframe with kinases as the index and datasets as the columns. The values in the dataframe should be either the activity (-log10(MW_p)) or fdr for each kinase in each dataset.
    result_type : str
        Type of hit to assess. 'rank' for activity, 'sig' for significance.
    k : int
        Rank cutoff for activity hits. Default is 10 (must be in top 10 most active). Only needed if result_type is 'rank'.
    p : float
        Significance cutoff for FPR hits. Default is 0.05. Only needed if result_type is 'sig'.
    """
    #if comparing multiple runs or to other algorithms, add to this list
    result_names = list(results.keys())

    ### where to store whether a specific kinase was a hit or missing
    hits_matrix = {'Y':pd.DataFrame(np.nan, index = range(meta['Y'].shape[0]), columns =['Kinase']+result_names),
               'ST':pd.DataFrame(np.nan, index = range(meta['ST'].shape[0]), columns =['Kinase']+result_names)}
    hits_matrix['Y']['Kinase'] = hits_matrix['Y']['Kinase'].astype(str)
    hits_matrix['ST']['Kinase'] = hits_matrix['ST']['Kinase'].astype(str)
    
    if result_type == 'sig':
        dataset_hits = {'Y':{'KSTAR':[],'KSEA':[],'PTM-SEA':[]},
               'ST':{}}
    
    ## where to store the total number of kinases with a prediction for each condition/algorithm
    num_kinases = {'Y':{},
                'ST':{}}
    dataset_hits = {'Y':{},
                'ST':{}}
    for name in result_names:
        num_kinases['Y'][name] = []
        num_kinases['ST'][name] = []
        dataset_hits['Y'][name] = []
        dataset_hits['ST'][name] = []


    #if comparing multiple KSTAR runs or to other 
    for alg in result_names:
        for mod in ['Y','ST']:
            test_number = 0
            #create variables for recording strictly accuracy for each modification type
            res = results[alg][mod]
            if result_type == 'rank':
                res =  results[alg][mod].rank(ascending = False)

            for i in meta[mod].index:
                #extact kinase and dataset name being tested from meta dataframe
                kin = meta[mod].loc[i, 'Kinases'].lstrip()
                dataset = meta[mod].loc[i, 'Dataset']

                #add kinase being tested to hits matrix
                hits_matrix[mod].loc[test_number, 'Kinase'] = kin
                
                #remove kinases without data (less relevant for KSTAR usually)
                tmp = res.dropna(subset = [dataset])

                #if perturbed kinase has predictions, check if algorithm correctly predicted kinase as differentially active
                num_kinases[mod][alg].append(tmp.shape[0])
                if kin in tmp.index:
                    if result_type == 'rank':
                        #check to see if score is in top10 and store whether it was a hit or miss
                        if res.loc[kin,dataset]<=k:
                            hits_matrix[mod].loc[test_number, alg] = 1
                        else:
                            hits_matrix[mod].loc[test_number, alg] = 0
                    elif result_type == 'sig':
                        test = results[alg][mod].loc[kin,dataset] <= p
                        if test:
                            hits_matrix[mod].loc[test_number, alg] = 1
                            dataset_hits[mod][alg].append(dataset+','+kin)
                        else:
                            hits_matrix[mod].loc[test_number, alg] = 0
                    else:
                        raise ValueError("result_type must be either 'rank' or 'sig'")
                                
                    test_number = test_number + 1

    if result_type == 'rank':
        return hits_matrix
    elif result_type == 'sig':
        return hits_matrix, dataset_hits
    
def calculate_Phit(hits_matrix):
    """
    Calculate Phit (average number of hits for each kinase across all datasets)

    Parameters
    ----------
    hit_matrix : pd.DataFrame
        Hit matrix from get_hit_matrix() function


    Returns
    -------
    grouped_hits : pd.DataFrame
        Averaged hits for each kinase
    """
    #calculate Phit for Y and ST
    Phit = pd.DataFrame(None, index = ['ST','Y'], columns = hits_matrix['Y'].columns)

    #Y
    #group by kinase and calculate the average Phit_k (fraction of times kinase is a hit when expected)
    hits_matrix_grouped = hits_matrix['Y'].groupby('Kinase').mean()

    #calculate Phit
    Phit.loc['Y'] = hits_matrix_grouped.sum()/hits_matrix_grouped.count()

    #ST
    #group by kinase and calculate the average Phit_k (fraction of times kinase is a hit when expected)
    hits_matrix_grouped = hits_matrix['ST'].groupby('Kinase').mean()

    #calculate Phit
    Phit.loc['ST'] = hits_matrix_grouped.sum()/hits_matrix_grouped.count()
    Phit = Phit.drop('Kinase', axis = 1)
    return Phit


def Phit_barplot(Phit_rank, Phit_sig, figsize = (6.5,3)):
    """
    Create barplot of average Phit for each implementation of KSTAR

    Parameters
    ----------
    Phit_rank : pd.DataFrame
        Phit values for KSTAR activity
    Phit_sig : pd.DataFrame
        Phit values for KSTAR FPR
    """
    fig, axes = plt.subplots(ncols = 2, figsize = figsize, sharey = 'row')
    fig.subplots_adjust(wspace = 0.05)

    #collapse Phit dataframes for plotting and combine
    Phit2 = Phit_sig.melt(ignore_index = False, value_name = '$P_{hit}$').reset_index()
    Phit2['Metric'] = 'Significance'
    Phit1 = Phit_rank.melt(ignore_index = False, value_name = '$P_{hit}$').reset_index()
    Phit1['Metric'] = 'Rank'
    plt_data = pd.concat([Phit1, Phit2], axis = 0).reset_index()
    plt_data.rename({'index':'Mod'}, axis = 1, inplace = True)
    colors = list(sns.color_palette('colorblind', n_colors = plt_data['variable'].nunique()))
    #plot tyrosine kinase results
    bar = sns.barplot(x = 'variable', y = '$P_{hit}$', hue = 'Metric', data = plt_data[plt_data['Mod'] == 'Y'], ax = axes[0])

    #adjust plot parameters
    axes[0].set_xlabel('')
    axes[0].set_ylabel('$\mathdefault{P_{hit}}$')
    axes[0].set_ylim([0,1])
    axes[0].get_legend().remove()
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation = 45, ha = 'right')
    axes[0].set_title('Tyrosine')
    #adjust colors of the bars so that they match the legend scheme
    alpha = [1,0.4]
    for j in range(len(alpha)):
        for i in range(len(colors)):
            if j == 0:
                bar.patches[i].set_color(colors[i])
                bar.patches[i].set_alpha(alpha[j])
            else:
                bar.patches[i+plt_data['variable'].nunique()].set_color(colors[i])
                bar.patches[i+plt_data['variable'].nunique()].set_alpha(alpha[j])

    #plot the serine/threonine barplot
    bar = sns.barplot(x = 'variable', y = '$P_{hit}$', hue = 'Metric', data = plt_data[plt_data['Mod'] == 'ST'], ax = axes[1])
    axes[1].set_xlabel('')
    axes[1].set_ylabel('')

    #legend
    legend_elements = [Patch(facecolor=colors[0], label = 'Rank ($\mathdefault{r\leq10}$)'),
                    Patch(facecolor=colors[0], alpha = 0.4, label='Significance ($\mathdefault{FDR\leq0.05}$)')]
    axes[1].legend(handles = legend_elements, bbox_to_anchor = (1,1))
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation = 45, ha = 'right')
    axes[1].set_title('Serine/Threonine')
    
    #adjust colors so that they match legend scheme
    alpha = [1,0.4]
    for j in range(len(alpha)):
        for i in range(len(colors)):
            if j == 0:
                bar.patches[i].set_color(colors[i])
                bar.patches[i].set_alpha(alpha[j])
            else:
                bar.patches[i+plt_data['variable'].nunique()].set_color(colors[i])
                bar.patches[i+plt_data['variable'].nunique()].set_alpha(alpha[j])
                

###########################################################################################################################################
######################## Robustness Analysis #############################################################################################
###########################################################################################################################################

def add_offcenter_grid(x = [0,5], y = [0,8], increment = 1, c = 'black', lw = 2):
    """
    Add bold grid around the triangles to delineate each experiment/algorithm

    Parameters
    ----------
    x : list
        x-axis range to plot grid
    y : list
        y-axis range to plot grid
    increment : int
        increment for grid lines
    c : str
        color of grid lines
    lw : int
        linewidth of grid lines

    """
    for i in range(x[0],x[1]+1):
        plt.plot([i,i],y, ls = '-', c = c, lw = lw)
    
    for j in range(y[0],y[1]+1):
        plt.plot(x, [j,j], ls = '-', c = c, lw = lw)

def kinase_specific_heatmap(hits_matrix_rank, hits_matrix_sig, mod = 'Y', figsize = (6,2)):
    """
    Plot heatmap to indicate Phit for individual kinases

    Parameters
    ----------
    hits_matrix_rank : pd.DataFrame
        Hit matrix for KSTAR activity rank
    hits_matrix_sig : pd.DataFrame
        Hit matrix for KSTAR activity significance
    mod : str
        Modification type to plot (Y or ST)
    figsize : tuple
        Size of figure to plot
    """
    #create matrix for average rank hits
    mat1 = hits_matrix_rank[mod].copy()
    mat1 = mat1.groupby('Kinase').mean()
    mat1.dropna(inplace = True)
    #create matrix for average significance hits
    mat2 = hits_matrix_sig[mod].copy()
    mat2 = mat2.groupby('Kinase').mean()
    mat2.dropna(inplace = True)
    mat2 = mat2.loc[mat1.index.values]

    plt.figure(figsize = figsize)
    N = mat1.shape[1]
    M = mat1.shape[0]
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)
    xs2,ys2 = np.meshgrid(np.arange(3,M+1), y)

    #create triangle objects to use for plotting
    triangles1 = [(i + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
    triangles2 = [(i + j*(M+1),i+1 + j*(M+1), i+1 + (j+1)*(M+1)) for j in range(N) for i in range(M)]
    triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1)
    triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2)

    #Plot upper left triangle (rank-based)
    img1 = plt.tripcolor(triang1, mat1.T.values.ravel(), cmap=plt.get_cmap('coolwarm'), vmax=1, vmin = 0, edgecolor = 'black', linewidths = 0.45)
    #Plot lower right triangle (significance-based)
    img2 = plt.tripcolor(triang2, mat2.T.values.ravel(), cmap=plt.get_cmap('coolwarm'), vmax=1, vmin = 0, edgecolor = 'black', lw = 0.45)

    #adjust ticks to match location and label of plot
    plt.xlim(x[0]-0.02, x[-1])
    plt.ylim(y[0], y[-1]+0.06)
    plt.yticks(np.arange(0.5, mat1.shape[1], 1),mat1.columns, rotation=0)
    plt.xticks(np.linspace(0.5,mat1.shape[0]-0.5,mat1.shape[0]), mat1.index)
    plt.xticks(rotation = 45, ha = 'right')
    #bold edges of each cell
    add_offcenter_grid(x = [0, mat1.shape[0]], y = [0,5])
    plt.title('Serine/Threonine')

def plot_top10_kinases(rankings, activities, cmap = sns.color_palette('colorblind'), kinase_to_highlight = None, study_labels = ['(8)', '(9)', '(10)', '(11)'], title = '',  figsize = (4, 9)):
    """
    For each dataset in the rankings dataframe from robustness analysis, plot the top 10 kinases by activity.

    Parameters
    ----------
    rankings : pd.DataFrame
        Rankings of kinases for each dataset in the experiment
    activities : pd.DataFrame
        Activity values for each kinase in each dataset
    cmap : list
        List of colors to use for plotting
    kinase_to_highlight : str or list
        Kinase(s) to highlight with unique color on the barplot (such as ABL1 for BCR-ABL fusion cell line, K562)
    study_labels : list
        List of labels to use for each dataset
    title : str
        Title to add to the plot
    figsize : tuple
        Size of the figure to plot
    
    """
    #set up figure
    fig, axes = plt.subplots(figsize = figsize, 
            nrows = rankings.shape[1]+1, ncols = 2,
            sharex = 'col',
            gridspec_kw = { 
                'height_ratios': [0.5]+[1 for i in range(rankings.shape[1])],
                'width_ratios': [1, 0.15]
            },)
    fig.subplots_adjust(wspace=0, hspace=0.1)

    for i in range(rankings.shape[1]):
        top10 = rankings.sort_values(by = rankings.columns[i])
        sample = top10.columns[i]
        kinases = top10.index[0:10]
        kinases = kinases[::-1]

        bar = axes[i+1,0].barh(kinases, -np.log10(activities.loc[kinases,sample]), color = cmap[0])
        plt.setp(axes[i+1,0].yaxis.get_majorticklabels(), fontsize = 10)

        #color kinases of interest on barplots
        if kinase_to_highlight is not None:
            if isinstance(kinase_to_highlight, str):
                loc = np.where(kinases == kinase_to_highlight)[0]
                if len(loc) > 0:
                    bar[loc[0]].set(color = cmap[3])
            elif isinstance(kinase_to_highlight, list):
                for kinase in kinase_to_highlight:
                    loc = np.where(kinases == kinase)[0]
                    if len(loc) > 0:
                        bar[loc[0]].set(color = cmap[3])

    

    axes[rankings.shape[1], 0].set_xlabel('Activity (-log10(p))')

    #Add title to first axes
    axes[0,0].annotate(title, (16.5,0.3), ha = 'center', fontsize = 16)
    axes[0,0].axis('off')

    #add study labels
    axes[0,1].axis('off')
    for i in range(rankings.shape[1]):
        axes[i+1, 1].annotate(study_labels[i], (0.5, 0.5), ha = 'center', va = 'center', fontsize = 14)
        axes[i+1,1].set_yticks([])
        #axes[i+1,1].axis('off')
    axes[rankings.shape[1],1].set_xticks([])

def robustness_dotplot(activities, fpr, context):
    #Setup subplots so that dendrograms are included
    fig, axes = plt.subplots(figsize = (10, 14), 
            nrows = 4, ncols = 2, 
            sharex = 'col', 
            sharey = 'row',
            gridspec_kw = {
                'height_ratios':[0.12,0.1,0.06, 1], 
                'width_ratios':[0.1,1]
            },)
    fig.subplots_adjust(wspace=0, hspace=0)



    dots = DotPlot(activities, 
                        fpr, 
                        legend_title='-log10(p-value)')
    #Cluster changes the sorting of the values array, so be sure to plot context last so that it is in the same sort.
    dots.drop_kinases_with_no_significance()
    dots.cluster(orientation = 'left', ax = axes[3,0], method='ward')
    dots.cluster(orientation = 'top', ax = axes[0,1], method='ward')
    dots.context(ax=axes[1,1],info = context, id_column = 'Sample', context_columns = ['Cell Line'], 
                orientation = 'top', margin = 300)

    axes[0,0].axis('off')
    axes[1,0].axis('off')
    axes[2,0].axis('off')
    axes[3,0].xaxis.set_visible(False)
    axes[2,1].margins(0.05,1)

    dots.dotplot(ax = axes[3,1])

    #create context barplot
    #labels = ['NSCLC' for i in range(7)] + ['CML' for i in range(4)]
    #studies = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11']
    #identify unique labels and assign colors
    #groups = np.unique(labels)
    #context_colors = sns.color_palette('colorblind', len(groups))
    #get xticks and bar length
    #ticks = axes[1,1].get_xticks()
    #len_bar = ticks[0]*2
    #color_order = list(range(len(groups)-1,-1,-1))
    #for lab, study,tick in zip(labels, studies, ticks):
    #    #decide color depending on label
    #    which_group = np.where(groups == lab)[0][0]
    #    axes[2,1].barh(0, len_bar, left = tick - len_bar/2, height = 0.25, color = context_colors[color_order[which_group]])
    #    axes[2,1].annotate(study, (tick, 0), ha = 'center', va = 'center', color = 'white')

def average_kinase_rankings(rankings):
    """
    Plot average rankings for NSCLC and CML datasets

    Parameters
    ----------
    rankings : pd.DataFrame
        Rankings of kinases for each dataset in the experiment
    """
    #plot average rankings
    fig, axes = plt.subplots(figsize = (3,3),
                        nrows = 1, ncols = 2,
                            sharey = 'row')
    fig.subplots_adjust(wspace = 0.1)

    #order rankings to group NSCLC and CML samples
    order = ['data:H3255_Guo', 'data:H3255_Rikova', 'data:HCC827_Guo', 'data:HCC827_Rikova', 'data:H3255_Moritz', 'data:HCC827_Beekhoff', 'data:H3255_Zhang', 'data:K562_Palma_2D', 'data:K562_Asmussen', 'data:K562_Beekhoff', 'data:K562_Palma_pYID']
    rankings = rankings[order]

    ##### get average ranks for nsclc samples #####
    #get nsclc sample results
    nsclc = rankings.iloc[:,0:7]
    averageRanks_nsclc = nsclc.mean(axis = 1)
    averageRanks_nsclc = averageRanks_nsclc.sort_values()


    ##### get average ranks for cml samples #####
    #get cml sample results
    cml = rankings.iloc[:,7:]
    averageRanks_cml = cml.mean(axis = 1)
    averageRanks_cml = averageRanks_cml.sort_values()

    colors = sns.color_palette('colorblind', 2)
    numKin = 5
    x = np.arange(numKin)/2
    ranks = averageRanks_nsclc
    high = 7
    bar = axes[0].bar(x,high - ranks[0:numKin].values, bottom = ranks[0:numKin],width = 0.4, color = colors[0])
    #plot first ranking for test case
    for i in range(len(x)):
        kinase = ranks.index[i]
        axes[0].annotate(kinase, (x[i], ranks.iloc[i]-0.25), ha = 'center', color = 'black',rotation = 90, fontsize = 11)

    plt.yticks(range(1,8))
    axes[0].set_ylabel('Average Activity Rank', fontsize = 12)
    #axes[0].set_xlim([-0.2,(numKin-1)/2+0.2])
    axes[0].set_xticklabels([])
    axes[0].set_ylim([0.1,high])

    ranks = averageRanks_cml
    high = 7
    bar = axes[1].bar(x,high - ranks[0:numKin].values, bottom = ranks[0:numKin],width = 0.4, color = colors[1])
    #plot first ranking for test case
    for i in range(len(x)):
        kinase = ranks.index[i]
        axes[1].annotate(kinase, (x[i], ranks.iloc[i]-0.25), ha = 'center', color = 'black',rotation = 90, fontsize = 11)

    axes[1].set_xticklabels([])
    axes[0].set_xlabel('NSCLC',fontsize = 12)
    axes[1].set_xlabel('CML', fontsize = 12)


    plt.gca().invert_yaxis()


##calculate jaccard indexes using the sites found in each datases
def jaccard(list1,list2):
    """
    Calculate jaccard similarity between two lists of items
    """
    intersection = len(list(set(list1) & set(list2)))
    union = len(list(set(list1) | set(list2)))
    jaccard = float(intersection)/union
    return jaccard

def phosphositeSimilarity(evidence, setsToCompare):
    """
    Calculate the similarity between two phosphoprotoemic datasets based on the sites identified in each experiment 
    (assessed using Jaccard similarity)
    
    Parameters
    ---------
    evidence: dataframe
        binary evidence indicating whether a site was identified in a given dataset (should contain 
        KSTAR_ACCESSION and KSTAR_SITE columns)
    setsToCompare: list of 2 strings
        indicates which columns (i.e. which experiments in evidence) to compare
    """
    evidence['descriptor'] = evidence['KSTAR_ACCESSION'] + '_' + evidence['KSTAR_SITE']
    #add 'data:' to front of columns so that it can be recognized
    setsToCompare = ['data:'+setsToCompare[i] for i in range(len(setsToCompare))]
    set1 = evidence[evidence[setsToCompare[0]] == 1]['descriptor'].values
    set2 = evidence[evidence[setsToCompare[1]] == 1]['descriptor'].values
    similarity = jaccard(set1,set2)
    return similarity
    
def makeTableComp_Jac(evidence, order = None):
    """
    Compare all similarlity between all experiments present in evidence (finds experiment columns based on 'data:' being
    present in column name)
    """
    #Use order given in evidence unless order parameter is given
    if order is not None:
        table = pd.DataFrame(None, columns = order, index = order)
    else:
        renamed_columns = [col.split(':')[1] for col in evidence.columns if 'data:' in col]
        table = pd.DataFrame(None, columns = renamed_columns, index = renamed_columns)
        
    for i in range(table.shape[0]-1):
        for j in range(i+1, table.shape[0]):
            set1 = table.columns[i]
            set2 = table.columns[j]
            setsToCompare = [set1,set2]
            similarity = phosphositeSimilarity(evidence,setsToCompare)
            #removeDataName from sets
            table.loc[set1,set2] = similarity
            table.loc[set2,set1] = similarity
            
    for i in range(table.shape[1]):
        table.iloc[i,i] = 1

    
    return table

def similarityHeatMap(similarity, labels = None, studies = None, significance = None, annot = None, legend_title = None, colorbar_label = 'Jaccard Similarity', context_colors = None, context_gap = 0, cmap = 'coolwarm', vmin = -1, vmax = 1,figsize = 10):
    """
    Plot triangle heatmap indicating the similarity between each study (either based on phosphosites identified or predicted 
    kinase activity).
    
    Parameters
    ----------
    similarity: pandas dataframe
        symmetric matrix indicating the similarity between all datasets (either jaccard similarity or spearmans rank)
    labels: list of length equal to number of experiments in similarity
        descriptive category to indicate coloring of bars (in this case, NSCLC vs. CML)
    studies: list of length equal to number of experiments in similarity
        indicates the name (or identifier) for each experiment to label on context barplots
    significance: pandas dataframe
        Provides significance (p-values) for each similarity value. If annot = 'Sig', must be included. Otherwise, does not affect
        plot
    annot: string
        Type of annotation to include on plot. 
            'Value': In each heatmap cell, also include the actual similarity value, in addition to color coding
            'Sig': In each heatmap cell, include both similarity value and signicance ('*' -> p<0.05)
            None: No annotation, just color coding
    """
    fig = plt.subplots(figsize = (figsize,figsize))
    size = (similarity.shape[0] + 2)*2
    #heatmap subplot
    ax1 = plt.subplot2grid((size,size+1), (0,2), rowspan = similarity.shape[0]*2+2, colspan = similarity.shape[1]*2+2)
    #colorbar subplot
    ax2 = plt.subplot2grid((size,size+1), (0,size), rowspan = similarity.shape[0]*2+2, colspan = 1)
    #context bars subplot
    ax3 = plt.subplot2grid((size,size+1), (0,0), rowspan = similarity.shape[0]*2+2, colspan = 2)
    #context bars subplot
    ax4 = plt.subplot2grid((size,size+1), (size-2,2), rowspan = 2, colspan = similarity.shape[0]*2+2)
    
    #Plot heatmap and colorbar
    #Create triangle matrix of booleans
    triangle = np.tril(np.ones(similarity.shape)).astype(bool)
    ##### 3 possible annotation options #######
    #1) In each cell indicate the value of similaritu
    if annot == 'Value':
        tri_df = similarity.where(triangle)
        heat = sns.heatmap(tri_df.astype(float), ax = ax1, cmap = cmap, vmin = vmin, vmax = vmax, annot = True, 
                           linewidths = 0.5, cbar_kws={'label': colorbar_label}, xticklabels = False, yticklabels = False, 
                           cbar_ax = ax2, annot_kws = {'size':7})
    # 2) In each cell indicate the value of similarity and whether it is significant (* -> p < 0.05)
    elif annot == 'Sig':
        #force to float if not already
        significance = significance.astype(float)
        vals = np.array([[str(round(similarity.iloc[i,j],2)) for j in range(similarity.shape[0])] for i in range(similarity.shape[1])])
        sig = np.array([[vals[i,j]+'*' if significance.iloc[i,j] < 0.05 else vals[i,j]+'' for j in range(significance.shape[1])] 
                        for i in range(significance.shape[0])])
        tri_df = similarity.where(triangle)
        heat = sns.heatmap(tri_df.astype(float), ax = ax1, cmap = cmap, vmin = vmin, vmax = vmax, annot = sig, fmt = '',
                           linewidths = 0.5, cbar_kws={'label': colorbar_label}, annot_kws = {'size':7},
                           xticklabels = False, yticklabels = False, cbar_ax = ax2)
        ax2.yaxis.label.set_size(10)
    # 3) Don't annotate
    else:
        tri_df = similarity.where(triangle)
        heat = sns.heatmap(tri_df.astype(float), ax = ax1, cmap = cmap, vmin = vmin, vmax = vmax,linewidths = 0.5, xticklabels = False, yticklabels = False, cbar_ax = ax2)
        
    #create context barplot
    #identify unique labels and assign colors
    groups = np.unique(labels)
    if context_colors is None:
        context_colors = sns.color_palette('colorblind', len(groups))

    top = True
    #determine the starting bottom position and the length of the bars needed
    bottom = similarity.shape[0]+1
    len_bar = bottom/similarity.shape[0]
    color_order = list(range(len(groups)-1,-1,-1))
    for lab, study in zip(labels, studies):
        #make first bar white/invisible since triangle barplot doesn't have anything on top
        bottom = bottom - len_bar
        #decide color depending on label
        which_group = np.where(groups == lab)[0][0]
        ax3.bar(0,len_bar-context_gap, bottom = bottom+context_gap, width = 0.2, color = context_colors[color_order[which_group]])
        ax3.annotate(study, (0,bottom+len_bar/2), ha = 'center', va = 'center', color = 'white',fontsize = 10)
            
    #bottom context bars need to be plotted in opposite order, but does similar as above
    top = True
    bottom = similarity.shape[0]+1
    for lab, study in zip(labels[::-1], studies[::-1]):
        #make first bar white/invisible since triangle barplot doesn't have anything on top
        bottom = bottom - len_bar
        #decide color depending on label
        which_group = np.where(groups == lab)[0][0]
        ax4.barh(0, len_bar-context_gap, left = bottom+context_gap, height = 0.3, color = context_colors[color_order[which_group]])
        ax4.annotate(study, (bottom+len_bar/2, 0), ha = 'center', va = 'center', color = 'white',fontsize = 10)

    
    # surround samples from the same study with boxes
    xy_list = [(0,2),(1,3),(5,9),(7,10)]
    for xy in xy_list:
        rect = patches.Rectangle(xy,0.98,0.98, fill = False, ec = 'black', lw = 1.5)
        ax1.add_artist(rect)

    #housekeeping
    ax1.axis('off')
    ax3.set_ylim([0,similarity.shape[0]+1])
    ax3.axis('off')
    #ax3.set_yticks([])
    #ax3.set_xticks([])
    
    ax4.set_xlim([0,similarity.shape[0]+1])
    ax4.axis('off')
    #ax4.set_yticks([])
    #ax4.set_xticks([])
    plt.axis('off')

def makeTableComp_Corr(evidence, order = None):
    """
    Calculate spearman rank correlations between kinase activity profiles obtained from each experiment (both r and p)
    """
    #find data columns
    if order is not None:
        evidence = evidence[order
                           ]
    data_columns = evidence.columns
    #Make table to store both the correlation and p-value
    table_r = pd.DataFrame(None, columns = data_columns, index = data_columns)
    table_p = pd.DataFrame(None, columns = data_columns, index = data_columns)
    for i in range(len(data_columns)-1):
        for j in range(i+1, len(data_columns)):
            set1 = data_columns[i]
            set2 = data_columns[j] 
            similarity, p = stats.spearmanr(evidence[set1],evidence[set2])
            table_r.loc[set1,set2] = similarity
            table_r.loc[set2,set1] = similarity
            table_p.loc[set1,set2] = p
            table_p.loc[set2,set1] = p
            
    for i in range(len(data_columns)):
        table_r.iloc[i,i] = 1
        table_p.iloc[i,i] = 1
    return table_r, table_p

def correlation_boxplot(melted, x_column = 'Method', y_column = 'Spearman Correlation', hue_column = 'Comparison Type', pal = sns.color_palette('colorblind'), figsize = (3,3)):
    """
    Given melted similarity dataframe with spearman correlation between activity profiles, plot boxplot with datapoints overlayed
    """
    plt.figure(figsize = figsize)
    pal = sns.color_palette('colorblind')
    with_col = pal[-1]
    bet_col = pal[4]
    palette = {'Within Tissues':with_col, 'Between Tissues':bet_col}
    swarm = sns.swarmplot(x= x_column, y = y_column, hue = hue_column, data=melted, dodge = True, palette = 'dark:black')
    box = sns.boxplot(x= x_column, y = y_column, hue = hue_column, data=melted, whis=np.inf, palette = palette)
    swarm.legend_.remove()

    #create legend
    patch1 = patches.Patch(color=with_col, label='Same Tissue Type')
    patch2 = patches.Patch(color=bet_col, label='Different Tissue Type')
    plt.legend(handles=[patch1,patch2], bbox_to_anchor = (0.95,1.25))


###########################################################################################################################################
######################## Sensitivity Analysis #############################################################################################
###########################################################################################################################################


def create_binary_evidence_for_sensitivity(test_data, greater = True, agg = 'mean', threshold = 1.0,  evidence_size = None):
    """
    Returns a binary evidence data frame from the sensitivity dataset (with subsets of original dataframe by targeted attack or random attack) according to the parameters passed in for method for aggregating
    duplicates and considering whether a site is included as evidence or not. Adapted from kinase activity class function, but applies different thresholds based on if the data is inhibition or

    Parameters
    ----------
    test_data: pd.DataFrame
        Dataframe containing the sensitivity data for the kinase activity
    greater: bool
        Indicates whether to take sites above or below the threshold (or higher or lower sites for evidence size). If true, will take the sites with higher quantification values. Default is True.
    agg: str
        Method for aggregating duplicates. Default is 'mean'.
    threshold: float
        Threshold for determining whether a site is included or not. Default is 1.0.
    evidence_size: int
        Number of sites to include in the evidence dataframe. If None, will include all sites above the threshold. Default is None.

    Returns
    -------
    evidence_binary : pd.DataFrame
        Matches the evidence dataframe of the kinact object, but with 0 or 1 if a site is included or not.
        This is uniquified and rows that are never used are removed.
    
        
    """
    #check to make sure datasets have data
    data_columns = [col for col in test_data.columns if 'data:' in col]


    #collapse sites into single row based on agg parameter
    test_data = test_data.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[data_columns].agg(agg).reset_index()
    
    #set the binary evidence for whether a site is included
    evidence_binary = test_data.copy()
    if evidence_size is None:
        for col in data_columns:
            if greater:
                evidence_binary[col] = evidence_binary[col] >= threshold
            else:
                evidence_binary[col] = evidence_binary[col] <= threshold
    else:
        for col in data_columns:
            #check how many non_nan sites there (if less than N, set n to be equal to number of sites available
            num_sites_available = evidence_binary.dropna().shape[0]
            if num_sites_available >= evidence_size:
                n = evidence_size
            else:
                n = num_sites_available
                print(f"{col} has less than {evidence_size} sites available, so using {n} sites instead")
                
            if greater:
                max_indices = np.argsort(-evidence_binary[col].values)[0:n]
                evidence_binary[col] = 0
                col_loc = np.where(evidence_binary.columns == col)[0][0]
                evidence_binary.iloc[max_indices, col_loc] = 1
            else:
                min_indices = np.argsort(evidence_binary[col].values)[0:n]
                evidence_binary[col] = 0
                col_loc = np.where(evidence_binary.columns == col)[0][0]
                evidence_binary.iloc[min_indices, col_loc] = 1

    #remove phosphorylation sites that were not selected in any experiment (useful for very large experiments where removing the need to copy data reduces time)
    evidence_binary.drop(evidence_binary[evidence_binary[data_columns].sum(axis=1) == 0].index, inplace = True) 

    #add back compendia/study bias information to binary evidence
    compendia = config.HUMAN_REF_COMPENDIA[['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS']]
    evidence_binary = evidence_binary.merge(compendia, on = ['KSTAR_ACCESSION', 'KSTAR_SITE'], how = 'left')
        
        
    return evidence_binary


import math
import random

def randomlySelect(df, fraction, testNum, dataset):
    """
    Randomly select a fraction of the sites (rows) from the dataframe. This is used to generate the datasets for the random attack tests.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to select from
    fraction : float
        The fraction of the rows from the dataframe to select
    testNum : int
        The number of the test, used to label the new column
    dataset : str
        The name of the dataset, used to label the new column
    
    Returns
    -------
    pandas.DataFrame
        A subset of the original dataframe containing the randomly selected rows
    """
    sample = df.sample(frac = fraction)
    sample.rename({f'data:All_{dataset}':f"data:random_{fraction}_test{testNum}_{dataset}"}, axis = 1, inplace = True)
    return sample

def studyBiasSelect(df, fraction, testNum,dataset):
    """
    Similar to randomlySelect, but selects sites based on the number of compendia they are found in, with the least studied sites selected first. To do so, we first shuffle the rows at random (introducing variability), sort by compendia, then grab the desired number of sites. This is used to generate the datasets for the targeted attack tests

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to select from
    fraction : float
        The fraction of the rows from the dataframe to select
    testNum : int
        The number of the test, used to label the new column
    dataset : str
        The name of the dataset, used to label the new column

    Returns
    -------
    pandas.DataFrame
        A subset of the original dataframe containing the study-biased selected rows

    """
    #shuffle order of dataframe, adding random element
    sample = df.sample(frac = 1)
    #sort based on num compendia
    sample = sample.sort_values(by = 'KSTAR_NUM_COMPENDIA', ascending = True)
    #pull sites equal to num_sites * fraction, starting with the least well studied sites
    numSites = round(df.shape[0]*fraction)
    sample = sample.iloc[0:numSites]
    
    #rename to mark the test
    sample.rename({f'data:All_{dataset}':f"data:biased_{fraction}_test{testNum}_{dataset}"}, axis = 1, inplace = True)
    return sample

def createTestDataset(df, data_col, fractions = [0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05], numTests = 5):
    """
    Create a dataset for testing sensitivity to data loss and study bias for a given column, considering different subsets of the original data. This function generates the datasets for the random and targeted attack tests.

    Parameters
    ----------
    df : pandas.DataFrame
        The mapped phosphoproteomic dataframe to select from
    data_col : str
        The column in the dataframe to test
    fractions : list of floats
        The fractions of the rows to select for the test datasets
    numTests : int
        The number of test datasets to generate for each fraction
    
    Returns
    -------
    pandas.DataFrame
        A dataframe containing the original data column, as well as the randomly and study-bias selected datasets for the given column. This data should be ready for KSTAR analysis
    """
    new_df = df[['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_PEPTIDE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS']+[data_col]].drop_duplicates(subset = ['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_PEPTIDE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS'])
    #rename data col to indicate 100%
    study = data_col.split(':')[1]
    new_df.rename({data_col:f'data:All_{study}'}, axis = 1, inplace = True)
    #save original version of dataframe
    orig_df = new_df.copy()
    for frac in fractions:
        for test in range(numTests):
            rand = randomlySelect(orig_df, frac, test, study)
            new_df = new_df.merge(rand, how = 'left')
            biased = studyBiasSelect(orig_df, frac, test, study)
            new_df = new_df.merge(biased, how = 'left')
            
    return new_df

def createMasterDatasets(data, fractions = [0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05], numTests = 5):
    """
    Create the master datasets for testing sensitivity to data loss and study bias for all columns in the dataframe. This function generates the datasets for the random and targeted attack tests.

    Parameters
    ----------
    data : pandas.DataFrame
        The mapped phosphoproteomic dataframe to select from

    Returns
    -------
    pandas.DataFrame
        A dataframe containing the original data columns, as well as the randomly and study-bias selected datasets for each column. This data should be ready for KSTAR analysis
    """
    data_cols = [col for col in data.columns if 'data:' in col]
    all_data = None
    for col in data_cols:
        col_data = data.dropna(subset = [col])
        if all_data is None:
            all_data = createTestDataset(col_data, col, fractions = fractions, numTests= numTests)
        else:
            tmp_data = createTestDataset(col_data,col, fractions = fractions, numTests = numTests)
            all_data = all_data.merge(tmp_data, how = 'outer', on = ['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_PEPTIDE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS'])
    
    return all_data

def get_sensitivity_data_to_test(test_data, meta, dataset_hits):
    """
    Given the sensitivity test data (generated by createMasterDatasets() function) and the datasets for which KSTAR correctly predicted the perturbed kinases (datasets_hits), extract which datasets KSTAR needs to be run, binning into stimulation and inhibition datasets.

    Parameters
    ----------
    test_data : pd.DataFrame
        The sensitivity test data generated by createMasterDatasets(), including all datasets in the benchmarking data
    meta : pd.DataFrame
        The metadata dataframe for the benchmarking data
    dataset_hits : pd.DataFrame
        The datasets for which KSTAR correctly predicted the perturbed kinases

    Returns
    -------
    test_data_up : pd.DataFrame
        The test data for the datasets where KSTAR correctly predicted the perturbed kinases for stimulation datasets
    test_data_down : pd.DataFrame
        The test data for the datasets where KSTAR correctly predicted the perturbed kinases for inhibition datasets
    """
    #need to split the dataset between up and down
    up_datasets = meta[meta['Direction'] == 'Up']
    down_datasets =meta[meta['Direction'] == 'Down']

    #### go through and only keep the datasets that KSTAR correctly predicted
    #### go through and only keep datasets where kstar succeeded, binning into stimulation and inhibition datasets #####
    #extract study name without 'data:'
    hits = np.unique([hit.split(',')[0].split(':')[1] for hit in dataset_hits])
    #iterate through hits and identify datasets to keep for stimulation and inhibition groups
    keep_up = []
    keep_down = []
    for col in test_data.columns:
        for h in hits:
            if h in col:
                for data in up_datasets:
                    if h in data:
                        keep_up.append(col)
                        break
                for data in down_datasets:
                    if h in data:
                        keep_down.append(col)
                        break
                break
    test_data_up = test_data[['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_PEPTIDE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS']+keep_up]
    test_data_down = test_data[['KSTAR_ACCESSION','KSTAR_SITE','KSTAR_PEPTIDE','KSTAR_NUM_COMPENDIA','KSTAR_NUM_COMPENDIA_CLASS']+keep_down]
    return test_data_up, test_data_down


def melt_sensitivity_results(fdr, dataset_hits, mod = 'Y', alg_name = 'KSTAR'):
    """
    Process the predictions on the sensitivity test data, extracting kinase predictions for the kinases of interest (those that are perturbed and were correctly identified by KSTAR with the full dataset), producing a long form dataframe with each replicate, fraction of data loss, and type of attack.

    Parameters
    ----------
    fdr : pd.DataFrame
        The predictions on the sensitivity test data
    dataset_hits : pd.DataFrame
        The datasets for which KSTAR correctly predicted the perturbed kinases
    mod : str
        The modification type to consider (either 'Y' or 'ST')
    alg_name : str
        The name of the algorithm used to make the predictions, used for labeling the rows of the dataframe in case comparisons are made across algorithms
    """
    # go through and extract information from tests
    res_all = None
    #iterate through each dataset hit and extract the results from the fpr data
    for hit in dataset_hits:
        #find columns associated with study, keep only those
        study = hit.split(',')[0].split(':')[1]
        keep = []
        for col in fdr.columns:
            if study in col:
                keep.append(col)
        res = fdr[keep]
        #extract only the kinase of interest
        kin = hit.split(',')[1]
        #melt the data into long format, and keep only kinases of interest
        res = res.melt(ignore_index = False).loc[kin]
        #add additional information about the dataset, kinase, fraction of data loss, type of selection
        res['Study'] = study
        res['Kinase'] = kin
        res['Fraction'] = res['variable'].apply(lambda x: 1 if 'All' in x else float(x.split('_')[1]))
        res['Selection'] = res['variable'].apply(lambda x: 'Biased' if 'biased' in x else ('random' if 'random' in x else 'All'))
        res.drop('variable', axis = 1, inplace = True)
        #add to any prior data
        if res_all is None:
            res_all = res.copy()
        else:
            res_all = pd.concat([res_all,res])

    #label with the algorithm and modification type (Y or ST) 
    res_all['Algorithm'] = alg_name
    res_all['Mod'] = mod
    return res_all

def calculate_tolerable_data_loss(results):
    """
    Calculate the tolerable data loss for the sensitivity test data, based on the fraction of replicates for which the FDR is less than or equal to 0.05. This is calculated for each study, kinase, fraction, and selection type. Initially store results in a nested dictionary, but reformat into a dataframe for plotting.

    Parameters
    ----------
    results : pd.DataFrame
        The results of the sensitivity test data, in long form (generated by melt_sensitivity_results() function). this should include both Y and ST data, and results for any algorithms or algorithm implementations you are testing

    Returns
    -------
    tdl_plt_data: pd.DataFrame
        A dataframe containing the tolerable data loss for each study, kinase, fraction, and selection type
    """
    #create nested dictionaries for storing the data
    tdl = {'random':{'Y':{},'ST':{}},
                                'Biased': {'Y':{},'ST':{}}}
    condition = {'random':{'Y':{},'ST':{}}, 'Biased': {'Y':{},'ST':{}}}

    tmp_all = results.copy()
    tmp_all['value'] = (tmp_all['value'] <= 0.05) * 1 #determine if FDR is <= 0.05, and prediction was a hit
    for mod in tdl['random'].keys():
        for alg in results['Algorithm'].unique():
            #add algorithm key to dictionaries
            if alg not in tdl['random'][mod].keys():
                tdl['random'][mod][alg] = []
                tdl['Biased'][mod][alg] = []
                condition['random'][mod][alg] = []
                condition['Biased'][mod][alg] = []

            #extract data for algorithm and modification type of interest
            tmp = tmp_all[tmp_all['Algorithm'] == alg]
            tmp = tmp[tmp['Mod'] == mod]
            tmp.dropna(subset = 'value', inplace = True)
            tmp = tmp.drop(['Algorithm','Mod'], axis = 1)

            #calculate the fraction of replicates for which FPR <= 0.05 for each study, kinase, fraction, and selection
            mean = tmp.groupby(['Study','Kinase','Fraction','Selection']).mean()
            mean = mean.reset_index()

            #for each study, kinase, and selection type, identify the point where FDR reaches above 0.05
            for study in mean['Study'].unique():
                study_data = mean[mean['Study'] == study]
                for kin in study_data['Kinase'].unique():
                    for selection in ['Biased','random']:
                        #grab data for specific kinase
                        condition[selection][mod][alg].append(study+','+kin)
                        test_data = study_data[study_data['Kinase'] == kin]
                        
                        #grab data specific to selection type and remove any nan
                        random = test_data[test_data['Selection'] != selection].copy()
                        random.dropna(subset = ['value'], inplace = True)

                        #starting with the largest subset, check if more than half of the replicates failed (i.e. FDR > 0.05). If they did, note this as the tolerable data loss point
                        for frac in [0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05]:
                            hits_frac = random.loc[random['Fraction'] == frac, 'value'].values
                            if hits_frac < 0.5:
                                tdl[selection][mod][alg].append(1-frac)
                                break
                            elif frac == random['Fraction'].min():
                                tdl[selection][mod][alg].append(1 - random['Fraction'].min())
                                break

    #reformat dictionary into dataframe
    tdl_plt_data = None
    for selection in ['Biased', 'random']:
        if selection == 'Biased':
            sel_name = 'Random'
        else:
            sel_name = 'Targeted'
        for mod in ['Y',"ST"]:
            for alg in tdl[selection][mod].keys():
                if tdl_plt_data is None:
                    tdl_plt_data = pd.DataFrame({'Tolerable Data Loss (%)':tdl[selection][mod][alg], 'Algorithm':alg, 'Mod': mod, 'Selection': sel_name, 'Condition': condition[selection][mod][alg]})
                else:
                    tdl_plt_data = pd.concat([tdl_plt_data, pd.DataFrame({'Tolerable Data Loss (%)':tdl[selection][mod][alg],'Algorithm':alg, 'Mod':mod, 'Selection': sel_name, 'Condition': condition[selection][mod][alg]})])
                    
    #convert to percentages
    tdl_plt_data['Tolerable Data Loss (%)'] = tdl_plt_data['Tolerable Data Loss (%)']*100

    return tdl_plt_data

def calculate_area_under_curve(curve_data, stop_fraction, x_col = 'Fraction', result_col = 'value'):
    """
    Calculate the area under the curve using the trapezoidal rule. Stop calculating the area when the fraction reaches the stop_fraction (no predictions are made beyond this point).
    """
    #make sure data is sorted from largest to smallest fraction of data
    curve_data = curve_data.sort_values(by = x_col, ascending = False)

    #calculate area under the curve with trapezoidal rule
    trap_sum = 0
    for frac in curve_data[x_col]:
        #check to see if passed stop fraction
        if frac < stop_fraction:
            break
        #continue calculating integral
        elif frac == 1 or frac == stop_fraction:
            trap_sum = curve_data.loc[curve_data[x_col] == frac, result_col].values + trap_sum
        else:
            trap_sum = 2*curve_data.loc[curve_data[x_col] == frac, result_col].values + trap_sum
    #trap area = (deltaX/2)* [f(1)+2*f(0.9)+...+f(xMin)]
    area = (5)/2*(trap_sum)
    return area

def calculate_sensitivity(results):
    """
    Calculate the sensitivity of the algorithm for each condition (study and kinase) based on the area under the curve of the biased results and the random results. The sensitivity to data loss is defined as the difference in the area under the curve between the random loss curve and the FDR obtained when using the full dataset. The sensitivity to study bias is defined as the difference between the area under the curve of the biased results and the area under the curve of the random results.

    Parameters:
    ----------
    results: pd.DataFrame
        The results of the sensitivity analysis, processed by melt_sensitivity_results() function. The dataframe should have the following columns: 'Study', 'Kinase', 'Fraction', 'Selection', 'Algorithm', 'Mod', 'value'. The 'value' column should contain the FDR values. The 'Selection' column should contain the selection type (random or biased). The 'Algorithm' column should contain the algorithm name. The 'Mod' column should contain the modification type (Y or ST). 
    """
    results2 = results.copy()
    results2['Condition'] =  results2['Study']+','+results2['Kinase']

    #missing values -> FPR = 1
    random_score = {'Y': {},'ST': {}}
    bias_score = {'Y': {},'ST': {}}
    diff_score = {'Y': {},'ST': {}}
    conditions = {'Y': {},'ST': {}}
    kinases = {'Y': {},'ST': {}}

    #focus only the curve from 0.95 to 0.5
    res2 = results2[results['Fraction'] >= 0.5]
    for mod in diff_score.keys():
        for alg in results['Algorithm'].unique():
            #add alg key to each dictionary
            if alg not in random_score[mod].keys():
                random_score[mod][alg] = []
                bias_score[mod][alg] = []
                diff_score[mod][alg] = []
                conditions[mod][alg] = []
                kinases[mod][alg] = []

            #extract relevant data (algorithm, and modificaiton type)
            tmp = res2[res2['Algorithm'] == alg]
            tmp = tmp[tmp['Mod'] == mod]
            tmp = tmp.drop(['Algorithm','Mod'], axis = 1)
            mean = tmp.groupby(['Study','Kinase','Fraction','Selection','Condition']).mean()
            mean = mean.reset_index()

            #iterate through each study and kinase, calculate the area under the curve for both random and biased results, as well as the difference
            for study in mean['Study'].unique():
                study_data = mean[mean['Study'] == study].copy()
                study_data.replace(np.nan, 1, inplace = True)
                for kin in study_data['Kinase'].unique():
                    #grab data of interest
                    test_data = study_data[study_data['Kinase'] == kin]
                    #get random and biased data: find the lowest fraction where both have predictions (will stop integral calculation here). Remove any np.nan
                    biased = test_data[test_data['Selection'] != 'random'].copy()
                    biased.dropna(subset = ['value'], inplace = True)
                    random = test_data[test_data['Selection'] != 'Biased'].copy()
                    random.dropna(subset = ['value'], inplace = True)
                    stop_fraction = np.max([random['Fraction'].min(), biased['Fraction'].min()])

                    #calculate the area under the curve of biased results via trapezoidal rule
                    biased_area = calculate_area_under_curve(biased, stop_fraction, x_col = 'Fraction', result_col = 'value')
                    bias_score[mod][alg].append(biased_area[0])

                    #calculate the area under the curve of random results via trapezoidal rule
                    random_area = calculate_area_under_curve(random, stop_fraction, x_col = 'Fraction', result_col = 'value')
                    base_area = 0.5*random.loc[random['Selection'] == 'All','value'].values
                    random_score[mod][alg].append(random_area[0])


                    diff_score[mod][alg].append(biased_area[0] - random_area[0])
                    conditions[mod][alg].append(test_data['Condition'].unique()[0])
                    kinases[mod][alg].append(test_data['Kinase'].unique()[0])

    #process data into plottable dataframe
    selection = 'Random'
    sensitivity = None
    mods = ['Y', 'ST']

    for mod in mods:
        for alg in random_score[mod].keys():
                if sensitivity is None:
                    sensitivity = pd.DataFrame({'Sensitivity':random_score[mod][alg], 'Algorithm':alg, 'Condition': conditions[mod][alg], 'Selection': 'Random', 'Mod':mod})
                else:
                    sensitivity = pd.concat([sensitivity, pd.DataFrame({'Sensitivity':random_score[mod][alg],'Algorithm':alg, 'Condition': conditions[mod][alg], 'Selection': 'Random','Mod':mod})])

    selection = 'Targeted'      
    for mod in mods:
        for alg in diff_score[mod].keys():
            sensitivity = pd.concat([sensitivity, pd.DataFrame({'Sensitivity':diff_score[mod][alg],'Algorithm':alg, 'Condition': conditions[mod][alg], 'Selection': 'Targeted', 'Mod':mod})])
                
    sensitivity.reset_index(inplace = True)

    return sensitivity

def plot_loss_curves(results, sensitivity, mod = 'Y', max_loss = 50, figsize = (6,6)):
    """
    For each kinase and algorithm, plot the average loss curve for the random and targeted attack results, with the sensitivity scores annotated on the plot.

    Parameters
    ----------
    results : pd.DataFrame
        The results of the sensitivity analysis, processed by melt_sensitivity_results() function. The dataframe should have the following columns: 'Study', 'Kinase', 'Fraction', 'Selection', 'Algorithm', 'Mod', 'value'. The 'value' column should contain the FDR values. The 'Selection' column should contain the selection type (random or biased). The 'Algorithm' column should contain the algorithm name. The 'Mod' column should contain the modification type (Y or ST).
    sensitivity : pd.DataFrame
        The sensitivity scores calculated by calculate_sensitivity() function. The dataframe should have the following columns: 'Sensitivity', 'Algorithm', 'Condition', 'Selection', 'Mod'.
    mod : str
        The modification type to consider (either 'Y' or 'ST')
    max_loss : float
        The maximum percent of data loss to consider. Default is 50.
    figsize : tuple
        The size of the figure to plot. Default is (6,6).
    """
    #split data into random and targeted attack results
    rdata = sensitivity[sensitivity['Selection'] == "Random"].copy()
    bdata = sensitivity[sensitivity['Selection'] == "Targeted"].copy()

    #extract kinase of interest for each condition
    rdata['Kinase'] = rdata['Condition'].apply(lambda x: x.split(',')[1])
    bdata['Kinase'] = bdata['Condition'].apply(lambda x: x.split(',')[1])

    #calculate the average sensitivity across all datasets associated with that kinase (if multiple)
    rand_sens = rdata.groupby(['Kinase', 'Algorithm'])['Sensitivity'].mean().reset_index()
    bias_sens = bdata.groupby(['Kinase','Algorithm'])['Sensitivity'].mean().reset_index()


    #grab results for modification of interest, and convert fraction to percentage data loss (inverse of fraction)
    tmp = results[results['Mod'] == mod].sort_values('Kinase')
    tmp['condition'] = tmp['Study']+','+tmp['Kinase']
    tmp['Fraction'] = (1-tmp['Fraction'])*100
    #restrict to max data loss
    tmp = tmp[tmp['Fraction'] <= max_loss]
    num_cases = tmp['Kinase'].nunique()

    #grab number of cases and the cases for each kinase and algorithm
    num_cases = tmp['Kinase'].nunique()
    cases = tmp['Kinase'].unique()
    algorithms = tmp['Algorithm'].unique()

    #plot
    fig, axes = plt.subplots(ncols = len(algorithms), nrows = num_cases, sharey = 'row', figsize = figsize)
    fig.subplots_adjust(wspace = 0.15)

    cases = tmp['Kinase'].unique()
    algorithms = results['Algorithm'].unique()
    for i in range(tmp['Kinase'].nunique()):
        kinase = cases[i]
        for j in range(len(algorithms)):
            tmp2 = tmp[tmp['Algorithm'] == algorithms[j]]
            tmp2 = tmp2[tmp2['Kinase'] == kinase]
            tmp2 = tmp2.drop(['Study', 'Kinase', 'Algorithm','Mod'], axis = 1)
            if tmp2.shape[0] > 0:
                mean = tmp2.groupby(['Fraction','Selection'])['value'].mean()
                mean = mean.reset_index()
                std = tmp2.groupby(['Fraction','Selection'])['value'].std()
                std = std.reset_index()
                bmean = mean[mean['Selection'] != 'random']
                bstd = std[std['Selection'] != 'random']
                rmean = mean[mean['Selection'] != 'Biased']
                rstd = std[std['Selection'] != 'Biased']
                if len(algorithms) == 1:
                    ax = axes[i]
                else:
                    ax = axes[i,j]

                ax.plot(bmean['Fraction'],bmean['value'], color = 'green', label = 'Biased')
                ax.plot(rmean['Fraction'],rmean['value'], label = 'Random')
                ax.fill_between(rmean['Fraction'], rmean['value'],rmean['value'].min(), color = 'blue', alpha = 0.2)

                #add annotations with sensitivity scores
                tmp3 = rand_sens[rand_sens['Algorithm'] == algorithms[j]]
                ax.annotate(round(tmp3.loc[tmp3['Kinase'] == kinase, 'Sensitivity'].values[0], 2), (0.98, 0.62), c = 'blue', fontsize = 9)
                #add annotations with sensitivity scores
                tmp4 = bias_sens[bias_sens['Algorithm'] == algorithms[j]]
                ax.annotate(round(tmp4.loc[tmp4['Kinase'] == kinase, 'Sensitivity'].values[0], 2), (0.98, 0.25), c = 'green', fontsize = 9)
                #axes[i,0].axhline(0.05, c = 'red', linestyle = 'dashed', alpha = 0.4)
                if bmean.dropna()['Fraction'].max() < max_loss:
                    coverage_loss = bmean.dropna()
                    x1 = coverage_loss['Fraction'].max()
                    barea_under_curve = bmean.replace(np.nan, 1)
                    new_row = {'Fraction':x1, 'Selection':'Biased','value':1}
                    barea_under_curve = barea_under_curve.append(new_row, ignore_index = True).sort_values('value', ascending = False).sort_values('Fraction', ascending = False)
                    rmean_extra = rmean.copy()
                    new_row = {'Fraction':x1, 'Selection':'Targeted','value':rmean.loc[rmean['Fraction'] == x1, 'value'].values[0]}
                    rmean_extra = rmean_extra.append(new_row, ignore_index = True).sort_values('value', ascending = False).sort_values('Fraction', ascending = False)
                    #add point where data is lost
                    ax.plot(x1, coverage_loss.loc[coverage_loss['Fraction'] == coverage_loss['Fraction'].max(), 'value'], marker = 'o', c = 'k', markersize = 3)
                    ax.plot([x1,x1], [coverage_loss.loc[coverage_loss['Fraction'] == x1, 'value'].values[0], 1], c = 'green', ls = 'dashed')
                    ax.fill_between(barea_under_curve['Fraction'], barea_under_curve['value'],rmean_extra['value'], color = 'green', alpha = 0.2)
                else:
                    ax.fill_between(bmean['Fraction'], bmean['value'],rmean['value'], color = 'green', alpha = 0.2)

            else:
                ax.annotate('N/A', (25, 0.4), ha = 'center', va = 'center', fontsize = 16)
            if algorithms[j] == 'KSTAR':
                ax.set_ylabel(kinase, rotation = 0, ha = 'right', va = 'center')
                
            ax.set_xlim([0,50])

            ax.set_ylim([-0.05,1])
            if i == 11:
                ax.set_xticks([0,10,20,30,40, 50])
                ax.set_yticks([])
            else:
                ax.set_xticks([])
                ax.set_yticks([])

    #add titles
    if len(algorithms) == 1:
        axes[0].set_title(algorithms[0])
    else:
        for i in range(len(algorithms)):
            axes[0,i].set_title(algorithms[i])

    #fig.supxlabel('Percent Data Loss')
    fig.supxlabel('Percent Data Loss (%)', y = 0.05)
    fig.supylabel('False Discovery Rate', x = -0.025)

def plot_tolerable_data_loss(tdl_plt_data, figsize = (3.75, 2.5)):
    """
    Plot the tolerable data loss obtained from calculate_tolerable_data_loss() function. The plot will show a separate point for each kinase and study.

    Parameters
    ----------
    tdl_plt_data : pd.DataFrame
        The tolerable data loss data, processed by calculate_tolerable_data_loss() function. The dataframe should have the following columns: 'Tolerable Data Loss (%)', 'Algorithm', 'Mod', 'Selection', 'Condition'.
    """
    #set up subplots
    fig, ax = plt.subplots(figsize = figsize, ncols = 2, sharey = 'row')
    #establish colors for plot
    blues = sns.color_palette('Blues', n_colors = 10)
    greens = sns.color_palette('Greens', n_colors = 10)
    palette = {'Random': 'blue', 'Targeted': 'green'}
    
    #get the number of algirthms that will plotted
    num_algs = tdl_plt_data['Algorithm'].nunique()
    algorithms = tdl_plt_data['Algorithm'].unique()
    plt_data_Y = tdl_plt_data[tdl_plt_data['Mod'] == 'Y']
    s1 = sns.swarmplot(x = 'Algorithm', y = 'Tolerable Data Loss (%)', hue = 'Selection', data = plt_data_Y, palette = palette, s= 3.5, dodge = True, alpha = 0.6, zorder = 0, ax = ax[0])
    fig.subplots_adjust(wspace = 0.1)
    ax[0].set_title('Tyrosine')
    ax[0].set_xlabel('')
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation = 35, ha = 'right')
    ax[0].set_ylim([0, 120])

    #plot serine threonine tdl
    plt_data_ST = tdl_plt_data[tdl_plt_data['Mod'] == 'ST']
    s2 = sns.swarmplot(x = 'Algorithm', y = 'Tolerable Data Loss (%)', hue = 'Selection',data = plt_data_ST, s = 3.5,palette = palette, alpha = 0.6, dodge = True,  ax = ax[1],zorder = 0)
    #adjust plot parameters
    ax[1].set_title('Serine/Threonine')
    ax[1].set_xlabel('')
    ax[1].set_ylabel('')
    ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation = 35, ha = 'right')
    ax[1].set_ylim([0, 120])

    #remove one legend
    s1.legend().remove()
    s2.legend(bbox_to_anchor = (1.05, 1))


    #add lines for medians (Y)
    for i in range(num_algs):
        #median for random attack
        median = plt_data_Y[(plt_data_Y['Selection'] == 'Random') & (plt_data_Y['Algorithm'] == algorithms[i])]['Tolerable Data Loss (%)'].median()
        ax[0].plot([i-0.4,i],[median,median], c = 'k', linewidth = 1.5)
        #median for targeted attack
        median = plt_data_Y[(plt_data_Y['Selection'] == 'Targeted') & (plt_data_Y['Algorithm'] == algorithms[i])]['Tolerable Data Loss (%)'].median()
        ax[0].plot([i,i+0.4],[median,median], c = 'k', linewidth = 1.5)

        #medians for ST
        #median for random attack
        median = plt_data_ST[(plt_data_ST['Selection'] == 'Random') & (plt_data_ST['Algorithm'] == algorithms[i])]['Tolerable Data Loss (%)'].median()
        ax[1].plot([i-0.4,i],[median,median], c = 'k', linewidth = 1.5)
        #median for targeted attack
        median = plt_data_ST[(plt_data_ST['Selection'] == 'Targeted') & (plt_data_ST['Algorithm'] == algorithms[i])]['Tolerable Data Loss (%)'].median()
        ax[1].plot([i,i+0.4],[median,median], c = 'k', linewidth = 1.5)

def plot_sensitivity(sensitivity, mod = 'Y', figsize = (1.88, 2.5), ylim1 = None, ylim2 = None):
    """
    Plot the sensitivity scores obtained from calculate_sensitivity() function. The plot will show a separate point for each kinase and study.

    Parameters
    ----------
    sensitivity : pd.DataFrame
        The sensitivity scores calculated by calculate_sensitivity() function. The dataframe should have the following columns: 'Sensitivity', 'Algorithm', 'Condition', 'Selection', 'Mod'.
    """
    plt_data_Y = sensitivity[sensitivity['Mod'] == mod]
    fig, ax = plt.subplots(ncols = 2, figsize = figsize)
    plt_data2 = plt_data_Y[plt_data_Y['Selection'] ==  'Random']

    algorithms = plt_data_Y['Algorithm'].unique()
    #get colors
    color = sns.color_palette('Blues', n_colors = 10)
    pal = {Algorithm: color[5] for Algorithm in plt_data_Y.Algorithm.unique()}
    #plot sensitivity to data loss
    plt_data2 = plt_data_Y[plt_data_Y['Selection'] ==  'Random']
    sns.boxplot(x = 'Selection', y = 'Sensitivity', hue = 'Algorithm', data = plt_data2,palette = pal, ax = ax[0])
    #adjust figure parameters
    ax[0].get_legend().remove()
    ax[0].tick_params(axis='y', which='both', colors='b')
    ax[0].tick_params(axis = 'x', rotation = 90)
    ax[0].set_ylabel('Sensitivity to Data Loss', c = 'blue')
    print(algorithms)
    if (len(algorithms) == 1):
        ax[0].set_xticklabels(algorithms)
    else:
        ax[0].set_xticks(np.arange(-0.5,0.5,len(algorithms)+1)[1:-1])
        ax[0].set_xticklabels(algorithms)

    ax[0].set_xlabel('')
    if ylim1 is not None:
        ax[0].set_ylim(ylim1)
    ax[0].set_xlim([-0.43,0.43])
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    #get colors
    color = sns.color_palette('Greens', n_colors = 10)
    pal = {Algorithm: color[5] for Algorithm in plt_data_Y.Algorithm.unique()}
    #plot sensitivity to study bias
    sns.boxplot(x = 'Selection', y = 'Sensitivity', hue = 'Algorithm', data = plt_data_Y[plt_data_Y['Selection'] ==  'Targeted'], palette = pal, ax = ax[1])
    #adjust plot parameters
    ax[1].set_ylabel('Sensitivity to Study Bias', rotation = 270, labelpad = 10, c = 'green')
    ax[1].get_legend().remove()
    ax[1].tick_params(axis='y', which='both', colors='g')
    ax[1].tick_params(axis = 'x', rotation = 90)
    if (len(algorithms) == 1):
        ax[1].set_xticklabels(algorithms)
    else:
        ax[1].set_xticks(np.arange(-0.5,0.5,len(algorithms)+1)[1:-1])
        ax[1].set_xticklabels(algorithms)
    ax[1].set_xlabel('')
    if ylim2 is not None:
        ax[1].set_ylim(ylim2)
    ax[1].set_xlim([-0.43,0.43])
    fig.subplots_adjust(wspace = 0)