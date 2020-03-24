import numpy as np
import bottleneck

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