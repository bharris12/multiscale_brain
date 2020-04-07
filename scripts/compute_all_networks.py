import pandas as pd
import numpy as np
import statsmodels.api as sm
import scanpy as sc
import bottleneck
from scipy import stats
import gc

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style='white', font_scale=1.25)
plt.rc("axes.spines", top=False, right=False)
plt.rc('xtick', bottom=True)
plt.rc('ytick', left=True)

from itertools import combinations
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)


import networkx as nx

import sys
sys.path.append('../scripts/')
sys.path.append('/home/bharris/Correlation_Coexpression/scripts/')
sys.path.append('/home/bharris/vshape/scripts/')

from rank import rank
from processify import processify

from argparse import ArgumentParser

def parse_args():

    parser = ArgumentParser(description='Decompose it!')

    parser.add_argument('--dataset',type=str,help='Which dataset to run on')

    return parser.parse_args()


def create_nw(data, replace_nans):
    """Compute Co-expression network from the data
    
    Core network building function. We always run with replace_nans = True
    Slicing single cell data will reguarly produce genes with no counts.
    And any correlation with a vector of all 0s is Nan.
    
    Arguments:
        data {np.array} -- Array of float values in shape of genes x cells
        replace_nans {bool} -- Flag for whether to replace Nans in network
    
    Returns:
        np.array -- ranked co-expression matrix of genes x genes 
    """
    nw = np.corrcoef(data)
    np.fill_diagonal(nw, 1)
    nw = rank(nw)
    if replace_nans:
        nw[np.isnan(nw)] = bottleneck.nanmean(nw)
    return nw


dataset_dict = pd.read_csv(
    '/home/bharris/biccn_paper/data/dataset_dict_biccn_sets_7.csv',
    index_col=0).to_dict()

genes = np.genfromtxt(
    '/home/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv',
    dtype=str)

args = parse_args()
dataset = args.dataset

logging.info(dataset)
andata = sc.read_h5ad(dataset_dict[dataset]['andata'])
sc.pp.normalize_total(andata, target_sum=1e6)
andata2 = andata[:, genes]
del andata
gc.collect()

andata = andata2[andata2.obs.class_label.isin(
    ['GABAergic', 'Glutamatergic'])]

del andata2
gc.collect()

expression = andata.to_df()
metadata = andata.obs

for column in ['class_label', 'subclass_label', 'cluster_label','joint_cluster_label']:
    logging.info(column)
    values = np.unique(metadata[column].values)

    for val in values:
        logging.info(val)
        mask = metadata[column] == val
        
        nw = create_nw(expression[mask].values.T, True)
        logging.info(f'{val} created')
        nw = pd.DataFrame(nw, index=genes, columns=genes)
        file_name = f'/home/bharris/biccn_paper/data/networks/{dataset}/coexpression_nw_{column}_{val}.hdf5'
        nw.to_hdf(file_name, 'nw')
        logging.info(file_name)
        del nw
        gc.collect()
del andata
gc.collect()