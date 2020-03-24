import mkl
mkl.set_num_threads(16)
import pandas as pd
import numpy as np

import sys
sys.path.append('/home/bharris/Correlation_Coexpression/scripts')
sys.path.append('/home/bharris/biccn_paper/scripts')
from egad import run_egad
import biccn_nw_perf_funcs as perf
from rank import rank

import gc
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

genes = perf.genes

def create_bootstrap_nw(datasets_selected):
    agg_nw = np.zeros([genes.shape[0], genes.shape[0]])
    for dataset in datasets_selected:
        logging.info(dataset)
        nw = pd.read_hdf(
            f'/home/bharris/biccn_paper/data/bulk_rna/networks/{dataset}_pearson_nw.hdf5',
            'nw')
        agg_nw += nw.values
        del nw
        gc.collect()
    return pd.DataFrame(rank(agg_nw), index=genes, columns=genes)

datasets_used = np.genfromtxt('/home/bharris/biccn_paper/data/bulk_rna/datasets_used.csv',dtype=str)
stephan_markers = pd.read_csv('/home/bharris/biccn_paper/data/stephan_marker_nw.csv',index_col=0)

n_bootstraps = 100
marker_bootstrap = []
for i in range(n_bootstraps):
    logging.info(i)
    nw = create_bootstrap_nw(np.random.choice(datasets_used, size=len(datasets_used)))
    egad_res = run_egad(stephan_markers, nw)[['AUC']]
    egad_res['iteration'] = i
    marker_bootstrap.append(egad_res)
    del egad_res, nw
    gc.collect()

marker_bootstrap = pd.concat(marker_bootstrap)
marker_bootstrap.to_csv('/home/bharris/biccn_paper/data/bulk_rna/bootstrap_100_markers.csv')