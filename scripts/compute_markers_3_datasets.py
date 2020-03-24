import pandas as pd
import numpy as np
from itertools import combinations

import sys
sys.path.append('/home/bharris/biccn_paper/scripts/')

import biccn_nw_perf_funcs as perf

dataset_dict = perf.dataset_dict
genes = perf.genes
all_datasets = np.array(list(dataset_dict.keys()))
nw_combinations = list(combinations(all_datasets, 4))
print(all_datasets)
print()
print()
for nw_ds in nw_combinations:
    marker_ds = all_datasets[~np.in1d(all_datasets, nw_ds)]
    ds_name = '_'.join(marker_ds)
    markers_nw = perf.make_marker_nw(perf.compute_markers(marker_ds), 100)
    print(ds_name)
    print()
    markers_nw.to_csv(
        f'/home/bharris/biccn_paper/data/markers_3ds/markers_{ds_name}.csv')
