import numpy as np
import pandas as pd
from scipy import stats, sparse
import bottleneck
import matplotlib.pyplot as plt

def run_egad_train_test(go, nw, **kwargs):
    """EGAD GeneList Train Test Functions:
    
    Wrapper to lower level functions for EGAD Train Test

    EGAD measures modularity of gene lists in co-expression networks. When running 
    the Train Test version, we measure how well one gene list predicts another.
    We compute every possible pairwise combination of gene pairs given in the GO

    
    The useful kwargs are: 
    int - {min,max}_count : limits for number of terms in each gene list, these are exclusive values


    Arguments:
        go {pd.DataFrame} -- dataframe of genes x terms of values [0,1], where 1 is included in gene lists
        nw {pd.DataFrame} -- dataframe of co-expression network, genes x genes
        **kwargs 
    
    Returns:
        pd.DataFrame -- dataframe of terms x terms where every value is an AUROC of training's (row) 
        term prediciton of the testing term (column)
    """
    assert nw.shape[0] == nw.shape[1] , 'Network is not square'
    assert np.all(nw.index == nw.columns) , 'Network index and columns are not in the same order'
    nw_mask = nw.isna().sum(axis=1) != nw.shape[1]
    nw = nw.loc[nw_mask, nw_mask].astype(float)
    np.fill_diagonal(nw.values, 1)
    return _runNV(go, nw, **kwargs)


def _runNV(go, nw, min_count=20, max_count=1000):

    #Make sure genes are same in go and nw
    genes_intersect = go.index.intersection(nw.index)

    go = go.loc[genes_intersect, :]
    nw = nw.loc[genes_intersect, genes_intersect]

    #Make sure there aren't duplicates
    duplicates = nw.index.duplicated(keep='first')
    nw = nw.loc[~duplicates, ~duplicates]

    go = go.loc[:, (go.sum(axis=0) > min_count) & (go.sum(axis=0) < max_count)]
    go = go.loc[~go.index.duplicated(keep='first'), :]

    roc = _new_egad(go.values, nw.values)

    #Put output in dataframe
    return pd.DataFrame(roc, index=go.columns, columns=go.columns)


def _new_egad(go, nw):

    nFold = go.shape[1]
    #Build Cross validated Positive
    x, y = np.where(go)

    ### NEW ####
    CVgo = {}
    for i in np.arange(nFold):
        fold = np.tile(go[:,i][:,None], nFold)
        CVgo[i] = fold
    CVgo = np.concatenate(list(CVgo.values()), axis=1)


    filtering = np.tile(go, nFold) #+ CVgo

    sumin = np.matmul(nw.T, CVgo)

    degree = np.sum(nw, axis=0)

    predicts = sumin / degree[:, None]

    np.place(predicts, CVgo > 0, np.nan)

    #Calculate ranks of positives
    rank_abs = lambda x: stats.rankdata(np.abs(x))
    predicts2 = np.apply_along_axis(rank_abs, 0, predicts)

    #Masking Nans that were ranked (how tiedrank works in matlab)
    predicts2[np.isnan(predicts)] = np.nan

    

    #negatives :filtering == 0
    #Sets Ranks of negatives to 0
    np.place(predicts2, filtering == 0, 0)

    #Sum of ranks for each prediction
    p = bottleneck.nansum(predicts2, axis=0)

    ##### NEW ######
   #Number of predictions
    #Number of 1's masked for each GO term for each CV
    n_p = np.sum(filtering, axis=0)

    #Number of negatives
    #Number of GO terms - number of postiive
    n_n = filtering.shape[0] - np.sum(filtering.astype(bool) | CVgo.astype(bool), axis=0)

    ###### NEW ######
    roc = (p / n_p - (n_p + 1) / 2) / n_n

    roc = roc.reshape(nFold, go.shape[1])
    
    
    return roc
