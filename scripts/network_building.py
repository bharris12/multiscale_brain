from rank import rank
import numpy as np
import pandas as pd
import bottleneck
import gc

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


def nw_aggregation(nw_paths, genes, file_key='nw'):
    """Function for aggregating co-expression networks
    
    Takes a list of paths to HDF5 files and reads in networks,
    avearges them and then re-ranks.

    Each HDF5 needs to be in the Pytable in the fixed format with
    the network stored under the key listed in the keyword argument file_key
    
    Arguments:
        nw_paths {list} -- list of strings or paths to HDF5 files
        genes {np.array} -- numpy array of genes for network
    
    Keyword Arguments:
        file_key {str} --  key in HDF5 network is stored under (default: {'nw'})
    
    Returns:
        pd.DataFrame -- Aggregate Network
    """

    agg_nw = np.zeros([genes.shape[0], genes.shape[0]])
    for nw_path in nw_paths:
        nw = pd.read_hdf(nw_path,file_key)
        fill = bottleneck.nanmean(nw.values,axis=None)
        agg_nw +=nw.loc[genes,genes].fillna(fill).values
        del nw
        gc.collect()

    return pd.DataFrame(rank(agg_nw),index=genes, columns=genes)