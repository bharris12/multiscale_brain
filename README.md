# Multiscale Co-expression in the Brain
Benjamin Harris, Megan Crow, Stephan Fischer, Jesse Gillis

Preprint : https://doi.org/10.1101/2020.03.31.018630

Publication: Under Review

 **These scripts, in their current form, will not be immediately usable**. We highlight the scripts that include critical functions for key algorithms and functions used in the analysis. We also include the notebooks used to create figures. 

### Environment
All of the analysis was done in either python3 or R. An overview of the python3 environment is in `python_requirements.txt`. In the repository  `r_session_info.yaml` describes the R environment. 

Unless otherwise stated, scripts are run on the command line and include mandatory command line arguments.

## Data 
Data to recapitulate these results can be aquired from [The NEMO archive](https://assets.nemoarchive.org/dat-ch1nqb7)

We used the aligned data and provided metadata for the 7 datasets. They were downloaded and stored as Andata objects (convert to Loom or SingleCellExperiments for use in R) along with the metadata. In the programs, datasets are referenced using a name and then the programs rely on a dictionary that is stored in a file to map the location of the Andata File and cluster annotation info. In the code this dictionary is stored as the variable `dataset_dict`. 

As github is not designed for sharing data, especially, larger data, we only include a few small and critical datafiles.

`datasets_used.csv` includes GEO numbers of the bulk expression data we used from GEMMA. The actual data can be downloaded and parsed using code in `download_gemma.r` and `parse_gemma_files.r` (Paths will need to be changed).

`highly_expressed_7_datasets_75k.csv` contains the list of 4,201 genes used in all of the analysis in the paper. We computed this list by computing the rank of the average expression of each gene in each dataset after excluding all non-neuronal cells. We then took the top 7,500 genes from each dataset and selected every gene that was in 6 of the 7 top 7,500 lists. The code for this can be found in `gene_selection_7_datasets.py`


## Results Figures 

In this repository there is a notebooks folder that includes ipython notebooks that were used to generate the results figures.
They are titled by the figures they include in them. The only figure not created in a notebook was the sankey (riverplot) in Figure 2c, the code for that figure is in the scripts folder in the file `figure2c_sankey.r`.


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
For immediate use the function is in `network_building.py`.

The rank function referenced in this script is in `rank.py` and ranks the values in the matrix between 0 and 1


We always use aggregate networks from each dataset and then aggregate datasets together. To build the dataset specific networks for the class, subclass and cluster level data we used `aggregate_networks.ipynb`. Below is the general network aggregation function, which is available in `network_building.py`

```
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
        agg_nw +=nw.values
        del nw
        gc.collect()

    return pd.DataFrame(rank(agg_nw),index=genes, columns=genes)
```

### Differential Expression

Meta-analytic marker lists are built by aggregating marker genes from individual datasets. For each dataset, a list of putative marker genes is built using a one-vs-all Wilcoxon test to prioritize genes (converted to equivalent AUROC score). The dataset-specific marker lists are then combined into meta-analytic lists ranked primarily by recurrence (number of datasets where gene was detected as DE), secondarily by auroc ("de") or fold change ("fc"). Note that subclass-specific markers are computed within classes, e.g. Vip markers are extracted by finding genes that are differentially expressed with respect to all other GABAergic subclasses.

To compute this the data was stored as (SingleCellExperiment)[https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html] Ojects. And run using 

After creating SingleCellExperimentObjects, changing paths in `datasets.R` run the below commands to compute the markers.
```
Rscript make_de_lists.R
Rscript make_meta_lists.R
```

For cases when less than 7 datasets are used for computing the markers we recompute reccurence for the subset of the datasets used.

### Performance of Modules 

To measure performance of Modules in the network we use the algorithm EGAD, originally published in [Ballouz et al 2017 Bioinformatics](https://doi.org/10.1093/bioinformatics/btw695).


![](https://github.com/bharris12/multiscale_brain/blob/master/figures/egad_cartoon.png)

We use two versions of EGAD, both exist in standalone files and can be imported into any workflow as is.
The files are 
```
egad.py
egad_train_test_terms.py
```
The standard EGAD algorithm is in the file `egad.py`, and is usable with any reference network (like GO) and any co-expression network. 

The function use it is `runNV`. Below is the doctring for that function

```def run_egad(go, nw, **kwargs):
    """EGAD running function
    
    Wrapper to lower level functions for EGAD

    EGAD measures modularity of gene lists in co-expression networks. 

    This was translated from the MATLAB version, which does tiled Cross Validation
    
    The useful kwargs are:
    int - nFold : Number of CV folds to do, default is 3, 
    int - {min,max}_count : limits for number of terms in each gene list, these are exclusive values


    Arguments:
        go {pd.DataFrame} -- dataframe of genes x terms of values [0,1], where 1 is included in gene lists
        nw {pd.DataFrame} -- dataframe of co-expression network, genes x genes
        **kwargs 
    
    Returns:
        pd.DataFrame -- dataframe of terms x metrics where the metrics are 
        ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    """
 ```
 
For figure 3g-h we use a variation of the original EGAD that computes pairwise EGAD scores between all gene lists. This version of EGAD is in the file `egad_train_test_terms.py` and is run using the function `run_egad_train_test`. The function does not compute anything meaningful for the diagonal, so for figure 3g we replace the diagonal with the normal EGAD output and in 3h we remove the diagonal.

Below is the docstring for the EGAD train test version.
```
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
```

The computation of the train_test egad creates a dense matrix of genes x terms^2, which can get very large in RAM for many terms. The most we ran it with was 200 terms. 

### Computing Metacells

Metacells were computed using the published [package](https://tanaylab.github.io/metacell/)

At the top of the `compute_metacells.sh` script lists the parameters used and the `metacell_scripts.r` file includes the exact functions and code used to run metacells. 

Other than the listed libraries in either the R session info file or at the top of the `metacell_scripts.r`, you just need a loom file for each dataset, that contains the class label column (to remove non-neuronal cells) and the `highly_expressed_7_datasets_75k.csv` file as the gene inputs. You also will have to change input and output paths to use. 

