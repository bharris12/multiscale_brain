import pandas as pd
import gc
import numpy as np
import sys

sys.path.append('/home/bharris/Correlation_Coexpression/scripts')

from rank import rank
from biccn_nw_perf_funcs import *


def load_individual_nws():

    nws_path = '/home/bharris/biccn_paper/data/networks/'

    class_label_pearson = {}
    subclass_label_pearson = {}
    cluster_label_pearson = {}
    compositional_pearson = {}
    joint_cluster_label_pearson = {}
    for dataset in dataset_dict:
        logging.info(dataset)
        nw = pd.DataFrame(rank(
            pd.read_hdf(
                f'{nws_path}{dataset}/pearson_agg_all_class_label.hdf5',
                'nw').values),
                          index=genes,
                          columns=genes)
        class_label_pearson[dataset] = nw
        del nw
        gc.collect()
        nw = pd.DataFrame(rank(
            pd.read_hdf(
                f'{nws_path}{dataset}/pearson_agg_all_subclass_label.hdf5',
                'nw').values),
                          index=genes,
                          columns=genes)
        subclass_label_pearson[dataset] = nw
        del nw
        gc.collect()
        nw = pd.DataFrame(rank(
            pd.read_hdf(
                f'{nws_path}{dataset}/pearson_agg_all_cluster_label.hdf5',
                'nw').values),
                          index=genes,
                          columns=genes)
        cluster_label_pearson[dataset] = nw
        del nw
        gc.collect()
        nw = pd.read_hdf(
            f'{nws_path}{dataset}/coexpression_nw_compositional.hdf5', 'nw')
        compositional_pearson[dataset] = nw
        del nw
        gc.collect()
        nw = pd.DataFrame(rank(
            pd.read_hdf(
                f'{nws_path}{dataset}/pearson_agg_all_joint_cluster_label.hdf5',
                'nw').values),
                          index=genes,
                          columns=genes)
        joint_cluster_label_pearson[dataset] = nw
        del nw
        gc.collect()

    pearson_nws = {
        'class_label': class_label_pearson,
        'subclass_label': subclass_label_pearson,
        'cluster_label': cluster_label_pearson,
        'joint_cluster_label': joint_cluster_label_pearson,
        'compositional': compositional_pearson
    }
    return pearson_nws


def load_all_networks():

    pearson_nws = load_individual_nws()
    all_datasets = np.array(list(dataset_dict.keys()))
    all_datasets_no_zeng_10x_nuc_v3 = np.delete(
        all_datasets,
        np.where(all_datasets == 'zeng_10x_nuc_v3')[0][0])

    #Compute Aggregate for all
    pearson_nws['class_label']['aggregate'] = combine_networks(
        all_datasets, 'class_label', pearson_nws)
    pearson_nws['subclass_label']['aggregate'] = combine_networks(
        all_datasets, 'subclass_label', pearson_nws)
    pearson_nws['cluster_label']['aggregate'] = combine_networks(
        all_datasets, 'cluster_label', pearson_nws)
    pearson_nws['compositional']['aggregate'] = combine_networks(
        all_datasets, 'compositional', pearson_nws)
    pearson_nws['joint_cluster_label']['aggregate'] = combine_networks(
        all_datasets, 'joint_cluster_label', pearson_nws)

    ##Compute Aggregaet no Zeng_10x_nuc_v3
    pearson_nws['class_label'][
        'aggregate_no_zeng_10x_nuc_v3'] = combine_networks(
            all_datasets_no_zeng_10x_nuc_v3, 'class_label', pearson_nws)
    pearson_nws['subclass_label'][
        'aggregate_no_zeng_10x_nuc_v3'] = combine_networks(
            all_datasets_no_zeng_10x_nuc_v3, 'subclass_label', pearson_nws)
    pearson_nws['cluster_label'][
        'aggregate_no_zeng_10x_nuc_v3'] = combine_networks(
            all_datasets_no_zeng_10x_nuc_v3, 'cluster_label', pearson_nws)
    pearson_nws['compositional'][
        'aggregate_no_zeng_10x_nuc_v3'] = combine_networks(
            all_datasets_no_zeng_10x_nuc_v3, 'compositional', pearson_nws)
    pearson_nws['joint_cluster_label']['aggregate_no_zeng_10x_nuc_v3'] = combine_networks(
        all_datasets_no_zeng_10x_nuc_v3, 'joint_cluster_label', pearson_nws)

    return pearson_nws
