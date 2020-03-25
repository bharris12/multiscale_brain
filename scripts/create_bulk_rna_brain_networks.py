import pandas as pd
import numpy as np
import statsmodels.api as sm
import scanpy as sc
import bottleneck
from scipy import stats
import gc
from sklearn.metrics import pairwise_distances

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



from argparse import ArgumentParser

def rank(data):
    """Rank normalize data
    
    Rank standardize data to make nonparametric
    
    Arguments:
        data {np.array} -- 2-D coexpression network
    
    Returns:
        np.array -- Rank normalized between 0 and 1 array
    """
    orig_shape = data.shape
    data = bottleneck.nanrankdata(data) - 1
    return (data / np.sum(~np.isnan(data))).reshape(orig_shape)

def parse_args():

    parser = ArgumentParser(description='Decompose it!')

    parser.add_argument('--dataset', type=str, help='Which dataset to run on')
    parser.add_argument('--metric',
                        type=str,
                        choices=['pearson', 'spearman', 'proportionality'],
                        required=True)
    parser.add_argument(
        '--genes-fn',
        default=
        '/home/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv',
        type=str,
        help='Name of file that cointains list of genes to use')
    parser.add_argument('--outname',
                        default='',
                        type=str,
                        help='Extra name to add to th output')
    parser.add_argument(
        '--outpath',
        default='/home/bharris/biccn_paper/data/bulk_rna/network/',
        help='Path to save networks in')
    return parser.parse_args()


def proportionality(x, y):
    num = bottleneck.nanvar(np.log1p(y) - np.log1p(x))

    denom = (bottleneck.nanstd(np.log1p(x)) +
             bottleneck.nanstd(np.log1p(y)))**2
    try:
        return num / denom
    except:
        return np.nan


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
    np.fill_diagonal(data, 1)
    nw = rank(nw)
    if replace_nans:
        nw[np.isnan(nw)] = bottleneck.nanmean(nw)
    return nw

args = parse_args()
genes = np.genfromtxt(args.genes_fn, dtype=str)

dataset = args.dataset
metric = args.metric

logging.info(dataset)
logging.info(metric)
dataset_path = f'/home/bharris/biccn_paper/data/bulk_rna/gemma_parsed_expression/{dataset}_expression.csv'

expression_data = pd.read_csv(dataset_path, index_col=0).loc[genes]
expression_data.index.name = 'genes'

#Deal with duplicate genes
expression_data = expression_data.reset_index().drop_duplicates(
    subset='genes').set_index('genes')
logging.info(expression_data.shape)
expression_data.dropna(thresh=int(genes.shape[0] / 2),
                       axis=1)  #Remove Samples with almost all Nans

if metric == 'pearson':
    nw = create_nw(expression_data.values, True)
elif metric == 'spearman':
    nw = create_nw(expression_data.rank().values, True)
elif metric == 'proportionality':
    expression_data = expression_data.dropna()  #Handle Nans in genes
    nw = pairwise_distances(expression_data.values,
                            metric=proportionality,
                            n_jobs=-1)

    #Put back genes
    nw = pd.DataFrame(nw,
                      index=expression_data.index,
                      columns=expression_data.index).loc[genes, genes].values
    nw = rank(nw)
    nw[np.isnan(nw)] = bottleneck.nanmean(nw)

nws_path = f'{args.outpath}{dataset}_{metric}_{args.outname}nw.hdf5'
pd.DataFrame(nw, index=genes, columns=genes).to_hdf(nws_path, 'nw')
