import json
import requests
import os
from pprint import pprint
from time import sleep
import pandas as pd

#below function obtained from KEA3 creators here: https://maayanlab.cloud/kea3/templates/api.jsp
def get_kea3_results(gene_set, query_name):
    """
    Function that takes a gene list, runs it through KEA3, and outputs results.
    Taken from KEA3 API at the site listed above.
    
    Parameters
    ----------
    gene_set: list
        contains all genes identified in an experiment
    query_name: string
        desired name of KEA3 run
        
    Returns
    -------
    dictionary containing all KEA3 results
    """
    ADDLIST_URL = 'https://amp.pharm.mssm.edu/kea3/api/enrich/'
    payload = {
        'gene_set': gene_set,
        'query_name': query_name
    }
    response = requests.post(ADDLIST_URL, data=json.dumps(payload))
    if not response.ok:
        raise Exception('Error analyzing gene list')
    sleep(1)
    return json.loads(response.text)

def uniprotToGene(acc):
    """
    Converts list of Uniprot accessions to gene names
    
    Parameters
    ----------
    acc: list
        list of Uniprot accessions for proteins/sites identified in experiment
        
    Returns
    -------
    list of gene names identified in experiment
    """
    import urllib.parse
    import urllib.request

    url = 'https://www.uniprot.org/uploadlists/'
    
    #create a query string for uniprot query
    queryString = ''
    for a in acc:
        queryString = "%s %s"%(queryString, a)
        
    params = {
    'from': 'ACC+ID',
    'to': 'GENENAME',
    'format': 'tab',
    'query': queryString
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
        
    uni_to_genes = {}

    for line in response.decode('utf-8').split('\n'):
        if line:
            l_arr = line.split('\t')
            uni_to_genes[l_arr[0]] = (l_arr[1])

    uni_to_genes.pop('From')


    #now walk through each row, create a unique index and add accession
    return list(uni_to_genes.values())


#rest are self created to operate with mapped KSTAR data
def runKEA3onDataset(df, data_cols = None, binarize = True, restart = False):
    """
    Runs KEA3 on all data columns in provided binarized (and mapped) 
    dataframe from a particular experiment
    
    Parameters
    ----------
    df: pandas dataframe
        binarized dataframe which indicates whether a particular protein/site are identified in experiment. Should be mapped and
        KSTAR_ACCESSION column.
    data_cols: list
        indicates which columns contain binarized evidence. If none, will look for columns with 'data:' in front
    """
    if data_cols is None:
        data_cols = [col for col in df.columns if 'data:' in col]
    
    binary = df.copy()
    if binarize:
        binary[data_cols] = (binary[data_cols] > 1)*1
    
    if os.path.isfile('temp_kea.tsv'):
        all_results = pd.read_csv('temp_kea.tsv', sep = '\t')
        data_cols = [col for col in data_cols if col not in all_results['Query Name'].unique()]
    else:
        all_results = None
        
    for col in data_cols:
        #get rows with 1 in data col
        tmp = binary[binary[col] == 1]
        #obtain kstar accessions for those cols
        acc = tmp['KSTAR_ACCESSION'].values
        #convert to gene name using uniprot query
        genes = uniprotToGene(acc)
        #feed through KEA3
        results = get_kea3_results(genes, col)
        #reformat results into dataframe with mean ranks
        results = pd.DataFrame(results['Integrated--meanRank'])
        #concat to previous results
        if all_results is None:
            all_results = results.copy()
        else:
            all_results = pd.concat([all_results,results], axis = 0)
            
        all_results.to_csv('temp_kea.tsv', sep = '\t')
        
    os.remove('temp_kea.tsv')
        
    return all_results
    
#rest are self created to operate with mapped KSTAR data
def runKEA3onDataset(df, data_cols = None, binarize = True, direction = 'Up'):
    """
    Runs KEA3 on all data columns in provided binarized (and mapped) 
    dataframe from a particular experiment
    
    Parameters
    ----------
    df: pandas dataframe
        binarized dataframe which indicates whether a particular protein/site are identified in experiment. Should be mapped and
        KSTAR_ACCESSION column.
    data_cols: list
        indicates which columns contain binarized evidence. If none, will look for columns with 'data:' in front
    """
    if data_cols is None:
        data_cols = [col for col in df.columns if 'data:' in col]
    
    binary = df.copy()
    if binarize:
        if direction == 'Up':
            binary[data_cols] = ((binary[data_cols] > 1)*1).values
        elif direction == 'Down':
            binary[data_cols] = ((binary[data_cols] < 1)*1).values
        else:
            print('needs to be either Up or Down for direction')
            exit()
    
    if os.path.isfile('temp_kea.tsv'):
        all_results = pd.read_csv('temp_kea.tsv', sep = '\t')
        data_cols = [col for col in data_cols if col not in all_results['Query Name'].unique()]
    else:
        all_results = None
        
    for col in data_cols:
        #get rows with 1 in data col
        tmp = binary[binary[col] == 1]
        #obtain kstar accessions for those cols
        acc = tmp['KSTAR_ACCESSION'].values
        #convert to gene name using uniprot query
        genes = uniprotToGene(acc)
        #feed through KEA3
        results = get_kea3_results(genes, col)
        #reformat results into dataframe with mean ranks
        results = pd.DataFrame(results['Integrated--meanRank'])
        #concat to previous results
        if all_results is None:
            all_results = results.copy()
        else:
            all_results = pd.concat([all_results,results], axis = 0)
            
        all_results.to_csv('temp_kea.tsv', sep = '\t')
    
    if os.path.isfile('temp_kea.tsv'):
        os.remove('temp_kea.tsv')
        
    return all_results