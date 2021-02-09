import numpy as np
import pandas as pd
import scanpy as sc
import gc

from egad import run_egad
from rank import rank

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

from argparse import ArgumentParser

##Functions

def load_and_preprocess_scData(fn):
    adata = sc.read(fn)
    adata = adata[adata.obs.class_label.isin(['GABAergic', 'Glutamatergic'])]
    sc.pp.normalize_total(adata, target_sum=1e6)
    adata = adata[:, genes]
    return adata
    
def compute_bulk_samples(expression_df, number_of_samples):
    bulk_samples = pd.DataFrame(index=adata.var_names)
    sample_size = int(expression_df.shape[0] / 20)
    i = 0
    while expression_df.shape[0] > sample_size:
        sample = expression_df.sample(n=sample_size)
        bulk_samples[i] = sample.sum()
        i += 1
        expression_df = expression_df.drop(sample.index)
        del sample
        gc.collect()
    return bulk_samples


def compute_bulk_nw(bulk_samples):
    nw = np.corrcoef(bulk_samples.values)
    rank(nw)
    return nw


def measure_performance(reference, nw, name, run_id):
    egad_res = run_egad(reference, nw)
    egad_res['reference'] = name
    egad_res['run_id'] = run_id
    return egad_res


##Args
run_id = 1
number_of_samples = 20
parser = ArgumentParser(description='Get Run ID')
parser.add_argument('--run-id',
                    type=int,
                    required=True,
                    help='Run ID for storing resutls')
parser.add_argument('--num-samples',
                    type=int,
                    default=20,
                    help='Number of samples to creates')
args = parser.parse_args()

number_of_samples = args.num_samples
run_id = args.run_id

##Load Data
dataset_dict = pd.read_csv(
    '/home/bharris/biccn_paper/data/dataset_dict_biccn_sets_7.csv',
    index_col=0).to_dict()
genes = np.genfromtxt(
    '/home/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv',
    dtype=str)

## Reference networks
go_mouse = pd.read_hdf('/home/bharris/GO_data/go_mouse_nw.hdf5', 'go')
go_slim = np.genfromtxt('/home/bharris/GO_data/aug4.GOslim', dtype=str)
go_mouse_slim = go_mouse[np.intersect1d(go_mouse.columns, go_slim)]

normal_markers = pd.read_csv(
    '/home/bharris/biccn_paper/data/stephan_marker_nw.csv', index_col=0)


##Compute Agg NW
agg_nw = np.zeros([genes.shape[0], genes.shape[0]])
datasets = np.array(list(dataset_dict.keys()))
np.random.shuffle(datasets)

for dataset in datasets:
    logging.info(dataset)

    adata = load_and_preprocess_scData(dataset_dict[dataset]['andata'])
    bulk_samples = compute_bulk_samples(adata.to_df(), number_of_samples)
    nw = compute_bulk_nw(bulk_samples)
    agg_nw += nw
    del adata, nw, bulk_samples
    gc.collect()

rank(agg_nw)
agg_nw = pd.DataFrame(agg_nw, index=genes, columns=genes)


## Compute Network Performance
logging.info("Running Egad on Aggregate Network")

egad_res_go = measure_performance(go_mouse_slim, agg_nw, 'GO', run_id)
egad_res_markers = measure_performance(normal_markers, agg_nw, 'Markers',
                                       run_id)
egad_res = pd.concat([egad_res_go, egad_res_markers])
out_fn = f'/home/bharris/biccn_paper/data/pseudobulk/pseudo_bulk_egad_res_{run_id}_n_samples_{number_of_samples}.csv.gz'
logging.info(out_fn)
egad_res.to_csv(
    out_fn)


