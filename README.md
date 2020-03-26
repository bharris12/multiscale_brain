# Multiscale Co-expression in the Brain
Scripts and some small datafiles for "Multiscale Co-expression in the Brain"

In order to do this analysis thousands of data files were created that are not feasible to share. 
There are lots of scripts and files needed to just manage all the paths and a large portion of all scripts are dedicated just to manage the files or parititoning data. 
As a result **these scripts, in their current form, will not be immediately usable**. 
However, in this README we highlight the scripts that include critical functions for key algorithms and functions used in the analysis. 

All of the analysis was done in either python3 or R. An overview of the python3 environment is in `python_requirements.txt`. In the repository  `r_session_info.yaml` describes the R environment. 

Most scripts are run on the command line and include command line arguments to change parameters. In general scripts are run on a per dataset basis. Datasets are referenced using a name and then the programs rely on a dictionary that is stored in a file to map the location of the Andata File and cluster annotation info. A few scripts, mainly related to the Network Performance section are used in the notebooks.

## Data 
Data to recapitulate these results can be aquired from [The NEMO archive](https://assets.nemoarchive.org/dat-ch1nqb7)

We used the aligned data and provided metadata for the 7 datasets. They were downloaded and stored as Andata objects (convert to Loom for use in R) along with the metadata. 

We do include some small critical datafiles here.

`datasets_used.csv` includes GEO numbers of the bulk expression data we used from GEMMA. The actual data can be downloaded and parsed using code in `download_gemma.r` and `parse_gemma_files.r` (Paths will need to be changed).

`highly_expressed_7_datasets_75k.csv` contains the list of 4,201 genes used in all of the analysis in the paper. We computed this list by computing the rank of the average expression of each gene in each dataset after excluding all non-neuronal cells. We then took the top 7,500 genes from each dataset and selected every gene that was in 6 of the 7 top 7,500 lists. The code for this can be found in `gene_selection_7_datasets.py`


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

Thousands of networks were built from all the different paritions of the data. We always use aggregate networks from each dataset and then aggregate datasets together. To build the dataset specific networks for the class, subclass and cluster level data we used `aggregate_networks.ipynb`. For aggregating the metacell networks we have that in the notebook `figure2bc_3bc.ipynb`. The script `compute_all_metacell_nws.py` is extremely RAM intensive for the larger datasets. We haven't be able to solve a memory leak issue, and it would likely be best to modify the script to create each network as it's own process and save to file, then read in by dataset all the intermediate networks. 

### Differential Expression

Differential expression to compute markers was done in R using the Mann Whitney test using custom scripts for efficicency. This is all computed on the entire list of genes in the data, not the 4,201 subset. To get a list of the top 100 markers for a given subclass we separate the data by class and then do 1 vs All DE for each subclass. The Mann Whiteny outputs an adjusted P-value, AUROC and log2FC. For each dataset we classify a gene as DE if the log2FC > 2 and the adjusted p value is <.05. For aggregating across datasets we take the average AUROC and compute recurrence, which is the number of datasets a gene is significantly DE. We rank the genes by recurrence then avg AUROC and take the top 100 genes within the 4,201 list. 

To compute this the data was stored as (SingleCellExperiment)[https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html] Ojects. And run using 

```
Rscript make_de_lists.R
Rscript make_meta_lists.R
```
To make work you will need to change the paths in `datasets.R` 

For cases when less than 7 datasets are used for computing the markers we recompute reccurence for the subset of the datasets used.

### Performance of Modules 

To measure performance of Modules in the network we use the algorithm EGAD, originally published in [Ballouz et al 2017 Bioinformatics](https://doi.org/10.1093/bioinformatics/btw695).


![](https://github.com/bharris12/multiscale_brain/blob/master/figures/egad_cartoon.png)

The standard EGAD algorithm is in the file `egad.py`.

The function to actually use it is `runNV`. Below is the doctring for that function

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
 
For figure 3g-h we use a special version of EGAD that computes pairwise EGAD scores between all gene lists. This version of EGAD is in the file `egad_train_test_terms.py` and is run using the function `run_egad_train_test`. The function does not compute anything meaningful for the diagonal, so for figure 3g we replace the diagonal with the normal EGAD output and in 3h we remove the diagonal.

The computation of the train_test egad creates a dense matrix of genes x terms^2, which can get very large in RAM for many terms. The most we ran it with was 200 terms. 

### Computing Metacells

Metacells were computed using the published [package](https://tanaylab.github.io/metacell/)

At the top of the `compute_metacells.sh` script lists the parameters used and the `metacell_scripts.r` file includes the exact functions and code used to run metacells. 

Other than the installed libraries, you just need a loom file for each dataset, that contains the class label column (to remove non-neuronal cells) and the `highly_expressed_7_datasets_75k.csv` file as the gene inputs. You also will have to change input and output paths to use. 

