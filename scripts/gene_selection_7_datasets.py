""" Script for calculating shared highly expressed genes 

Script calculates top 7,500 highly expressed genes for each dataset and takes reccurent 
genes across the 7 datasets

It uses genes that are reccurent in 6 of the 7 datsets. 
And also only computes it using good neuron cells

The 7 datsets are:
zeng_10x_cell_v2
zeng_10x_nuc_v2
zeng_smart_cell
zeng_smart_nuc
zeng_10x_cell_v3
zeng_10x_nuc_v3
macoscko_10x_nuc_v3

"""
import pandas as pd
import numpy as np
import gc
import sys
sys.path.append('/home/bharris/vshape/scripts/')

import gene_selection_replicability as hvg

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

dataset_dict = pd.read_csv('/home/bharris/biccn_paper/data/dataset_dict_biccn_sets_7.csv',
                           index_col=0).to_dict()

selected_genes = []
for dataset in dataset_dict:
    logging.info(dataset)
    andata = hvg.read_and_process_andata(dataset_dict[dataset]['andata'],
                                         obs_col='class_label',
                                         obs_value1='GABAergic',
                                         obs_value2='Glutamatergic')
    hvgs = hvg.select_genes(andata, 'highly_expressed', n_top_genes=7500)
    hvgs.name = dataset
    selected_genes.append(hvgs)

gene_matrix = pd.concat(selected_genes, axis=1)
gene_matrix.fillna(0, inplace=True)

mostly_reccurent = gene_matrix.sum(axis=1) > len(dataset_dict) - 1

logging.info(mostly_reccurent.sum())
np.savetxt('/home/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv',
           mostly_reccurent.index.values[mostly_reccurent.values],
           fmt='%s')
