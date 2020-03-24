import pandas as pd
import numpy as np
import logging
from itertools import combinations
import re
import gc
from scipy import stats
from scipy import io
import bottleneck

import sys

sys.path.append('/home/bharris/Correlation_Coexpression/scripts/')
sys.path.append('/home/bharris/vshape/scripts/')

from rank import rank
from egad import run_egad

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

dataset_dict = pd.read_csv(
    '/home/bharris/biccn_paper/data/dataset_dict_biccn_sets_7.csv',
    index_col=0).to_dict()
genes = np.genfromtxt(
    '/home/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv',
    dtype=str)
subclass_recurrence = pd.read_csv(
    '/home/bharris/biccn_paper/data/subclass_recurrence.csv', index_col=0)

def evaluate_all_networks(gene_lists, metric_nws):
    results = []
    for metric in metric_nws:
        logging.info(metric)
        for nw_name in metric_nws[metric]:
            logging.info(nw_name)
            res = run_egad(gene_lists,
                           metric_nws[metric][nw_name],
                           min_count=10)[['AUC']]
            res['Metric'] = metric
            res['nw_name'] = nw_name
            results.append(res)
    return results

def combine_networks(nw_list, metric, metric_nws):
    agg = np.zeros([genes.shape[0], genes.shape[0]])
    for nw in nw_list:
        agg += metric_nws[metric][nw].values
    return pd.DataFrame(rank(agg), index=genes, columns=genes)

def fix_ds_list(ds_list):
    #Stephan's file names are slightly different than mine
    if 'zeng_10x_cell' in ds_list:  #Switching to Stephan's Annotations
        ds_list[ds_list == 'zeng_10x_cell'] = 'zeng_10x_cell_v2'
    if 'zeng_10x_nuc' in ds_list:
        ds_list[ds_list == 'zeng_10x_nuc'] = 'zeng_10x_nuc_v2'
    substitue_nuc = np.vectorize(lambda x: re.sub('nuc', 'nuclei', x))
    substitute_cell = np.vectorize(lambda x: re.sub('cell', 'cells', x))

    ds_list = substitute_cell(
        substitue_nuc(ds_list))  #Stephan writes out nuceli
    return ds_list

def compute_recurrent(marker_stats, dataset):
    #Compute markers with
    recur = ((marker_stats.fdr < .05).astype(float) *
             (marker_stats.fold_change > 2))
    recur = recur.reset_index().drop_duplicates(
        subset='gene').set_index('gene')[0][genes]
    recur.name = dataset
    return recur

def compute_markers(ds_list):
    selected_subclasses = subclass_recurrence.index[
        subclass_recurrence[ds_list].sum(axis=1) == ds_list.shape[0]].values
    stephan_path = '/tyronedata/fischer/de_lists/'

    ds_list = fix_ds_list(ds_list)
    substitue_underscore = np.vectorize(lambda x: re.sub('_', '.', x))

    #Also compute GABA and Glut for everyone
    selected_subclasses = np.append(selected_subclasses,
                                    ['GABAergic', 'Glutamatergic','Non.Neuronal'])
    gene_ranks = []
    for subclass in selected_subclasses:
        logging.info(subclass)
        recurrence = []
        auroc = []
        for dataset in ds_list:
            try:
                marker_stats = pd.read_csv(
                    f'{stephan_path}dataset_specific/{dataset}_{substitue_underscore(subclass)}.txt',
                    comment='#',
                    index_col=0,
                    delimiter=' ')
                recur = compute_recurrent(marker_stats, dataset)
                recurrence.append(recur)
                roc = marker_stats.loc[genes, 'auroc']
                roc.name = dataset
                auroc.append(roc)
            except:
                logging.info(f'{subclass} not in {dataset}')
        if len(auroc) < ds_list.shape[0]:
            logging.info(f'Skipping {subclass}')
            continue
        auroc = pd.concat(auroc, axis=1).mean(axis=1)
        recurrence = pd.concat(recurrence, axis=1).sum(axis=1)
        genes_r = pd.Series(pd.concat([auroc, recurrence], axis=1).sort_values(
            [1, 0], ascending=False).index.values,
                            name=subclass)
        gene_ranks.append(genes_r)
    return pd.concat(gene_ranks, axis=1)



def make_marker_nw(markers, n_genes):
    nw = pd.DataFrame(0, index=genes, columns=markers.columns)
    for mark in markers.columns:
        nw.at[markers[mark].values[:n_genes], mark] = 1
    return nw

def add_attributes(df, metric, name):
    df['metric'] = metric
    df['name'] = name
    df.index.name = 'labels'
    return df