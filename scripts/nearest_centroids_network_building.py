import mkl
mkl.set_num_threads(16)

import pandas as pd
import numpy as np
import scanpy as sc
import bottleneck
from scipy import stats, io, spatial

import gc
import re

import sys
sys.path.append('../scripts/')
sys.path.append('/home/bharris/Correlation_Coexpression/scripts/')
sys.path.append('/home/bharris/vshape/scripts/')

from rank import rank
from egad import run_egad
import biccn_nw_perf_funcs as perf

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

from argparse import ArgumentParser

genes = perf.genes
dataset_dict = perf.dataset_dict


def parse_args():
    parser = ArgumentParser(description='Decompose it!')
    parser.add_argument('--dataset', type=str, help='Name of dataset to test')
    parser.add_argument(
        '--level',
        choices=['class', 'subclass', 'cluster', 'joint-cluster'],
        required=True,
        help='Level of Classification to test')
    parser.add_argument('--n_cells',
                        default=50,
                        help='Number of cells to use in network')
    return parser.parse_args()


def read_dataset(dataset, level, level_values):
    andata = sc.read_h5ad(dataset_dict[dataset]['andata'])
    sc.pp.normalize_total(andata, target_sum=1e6)
    andata = andata[:, genes]

    andata = andata[andata.obs[level].isin(level_values)]
    return andata


def select_cells(andata, level, n_cells, centroids):
    included_barcodes = []
    for a, b in andata.to_df().groupby(andata.obs[level]).groups.items():
        included_barcodes.extend(b[stats.rankdata(
            spatial.distance.cdist(andata.to_df().loc[b].values, centroids[
                [a]].values.T)) <= n_cells])
    andata = andata[included_barcodes]
    return andata


def create_nw(data, replace_nans=True):
    array = data.values.T
    nw = np.corrcoef(array)
    np.fill_diagonal(array, 1)
    nw = rank(nw)
    if replace_nans:
        nw[np.isnan(nw)] = bottleneck.nanmean(nw)
    return nw


def main():
    args = parse_args()
    dataset = args.dataset
    if args.level == 'class':
        level = 'class_label'
        level_values = ['GABAergic', 'Glutamatergic']
    elif args.level == 'subclass':
        level = 'subclass_label'
        spec = pd.read_csv(
            '/home/bharris/biccn_paper/data/networks/nearest_centroid/subclass_spec_table.csv',
            index_col=0)
        level_values = spec.columns[spec.loc[args.dataset].astype(bool)].values
    elif args.level == 'joint-cluster':
        level = 'joint_cluster_label'
        spec = pd.read_csv(
            '/home/bharris/biccn_paper/data/networks/nearest_centroid/joint_cluster_spec_table.csv',
            index_col=0)
        level_values = spec.columns[spec.loc[args.dataset].astype(bool)].values
    else:
        level = 'cluster_label'
        spec = pd.read_csv(
            '/home/bharris/biccn_paper/data/networks/nearest_centroid/cluster_spec_table.csv',
            index_col=0)
        level_values = spec.columns[spec.loc[args.dataset].astype(bool)].values

    logging.info(f'{dataset} : {level} ')
    andata = read_dataset(args.dataset, level, level_values)
    logging.info(f'{dataset} : {level} : Computing Centroids')
    centroids = andata.to_df().groupby(andata.obs[level]).mean().T

    logging.info(f'{dataset} : {level} : Selecting Cells')
    andata = select_cells(andata, level, args.n_cells, centroids)

    logging.info(f'{dataset} : {level} : Computing Networks')
    networks = andata.to_df().groupby(andata.obs[level]).apply(create_nw)
    path = '/home/bharris/biccn_paper/data/networks/nearest_centroid/'
    for name in networks.index:
        logging.info(f'{dataset} : {level} : writing {name}')
        pd.DataFrame(networks[name], index=genes, columns=genes).to_hdf(
            f'{path}{args.dataset}/nearest_centroid_nw_{level}_{name}_{args.n_cells}.hdf5',
            'nw')


main()
