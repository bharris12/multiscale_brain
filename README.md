# Multiscale Co-expression in the Brain
Scripts and some small datafiles for "Multiscale Co-expression in the Brain"

In order to do this analysis thousands of data files were created that are not feasible to share. 
There are lots of scripts and files needed to just manage all the paths and a large portion of all scripts are dedicated just to manage the files or parititoning data. 
As a result most of these scripts, in their current form will not be directly usable. 
However, in this README we highlight the scripts that include functions for key algorithms and functions used in the analysis. 

## Data 
Data to recapitulate these results can be aquired from [The NEMO archive](https://assets.nemoarchive.org/dat-ch1nqb7)

We used the aligned data and provided metadata for the 7 datasets. They were downloaded and stored as Andata objects (convert to Loom for use in R) along with the metadata. 


## Results Figures 

In this repository there is a notebooks folder that includes ipython notebooks that were used to generate the results figures.
They are titled by the figures they include in them. The only figure not created in a notebook was the sankey (riverplot) in Figure 2c,
the code for that figure is in the scripts folder in the file `figure2c_sankey.r`.


## Algorithms
In our analysis there are 4 main algorithms/functions used in the analysis 

1. Network Building
2. Differential Expression
3. Performance of modules within the network (EGAD)
4. Computing Metacells

We will highlight the specific code necessary for each algoirthm/function listed above

### Network Building
The core network building function is 
```
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
    np.fill_diagonal(data, 1)
    nw = rank(nw)
    if replace_nans:
        nw[np.isnan(nw)] = bottleneck.nanmean(nw)
    return nw
 ```
This function was used in `compute_all_networks.py`, `compute_all_metacell_nws.py` and `create_bulk_rna_brain_nws.py`

The rank function referenced in this script is in `rank.py` and ranks the values in the matrix between 0 and 1


### Differential Expression

### Performance of Modules 

### Computing Metacells



