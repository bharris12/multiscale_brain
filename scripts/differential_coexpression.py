import numpy as np
import gc
dropna = lambda x: x[np.isfinite(x)]

def build_null(n_datasets, n_genes):
    return np.random.uniform(
        size=[n_genes, n_genes, n_datasets]).mean(axis=2)

def calculate_quantile_fdr(empirical,
                           null,
                           n_bins=1000000,
                           center=False,
                           scaling_factor=1):
    """Compute FDR by putting each gene pair into tiny bins and computing for quanties
    
    Arguments:
        empirical {[type]} -- [description]
        null {[type]} -- [description]
     
    Keyword Arguments:
        n_bins {number} -- Number of Bins, more is better (default: {1000000})
        center {bool} --  Whether to center and flip networks to make two sided (default: {False})
        scaling_factor {number} -- If the dimensions of the null are greater than the dimensions of the empirical (default: {1})
    
    Returns:
        [type] -- [description]
    """
    if empirical.ndim == 2:
        np.fill_diagonal(empirical, np.nan)
        np.fill_diagonal(null, np.nan)

    #Center for two sided
    if center:
        center_empirical = np.abs(empirical - .5)
        center_null = np.abs(null - .5)
    else:
        center_empirical = empirical
        center_null = null

    null_quantiles = np.nanquantile(center_null, np.linspace(0, 1, n_bins))

    empirical_hist = np.histogram(dropna(center_empirical), null_quantiles)
    null_hist = np.histogram(dropna(center_null), null_quantiles)

    fdrs = np.flip(
        np.cumsum(np.flip(null_hist[0]) / scaling_factor) / np.cumsum(
            np.flip(empirical_hist[0])))
    

    empirical_digit = np.digitize(dropna(center_empirical), null_quantiles)

    empirical_digit[empirical_digit ==
                    empirical_digit.max()] = empirical_digit.max() - 1

    quantile_fdr = np.ones_like(center_empirical)

    quantile_fdr[np.isfinite(center_empirical)] = fdrs[empirical_digit - 1]
    quantile_fdr[np.isnan(center_empirical)] = np.nan

    empirical_density = np.flip(np.cumsum(np.flip(
        empirical_hist[0]))) / np.isfinite(center_empirical).sum()

    quantile_fdr[quantile_fdr > 1] = 1
    fdrs[fdrs>1]  =1 
    del empirical_hist, center_empirical, null_quantiles, center_null, empirical_digit
    gc.collect()

    return quantile_fdr, fdrs, empirical_density
