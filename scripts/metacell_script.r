library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(metacell)
library(loomR)
library(getopt)

rm(list = ls())

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

spec = matrix(
  c(
    'loomfile'     ,
    'l',
    1,
    "character",
    'genefile'     ,
    'g',
    1,
    "character",
    'dataset_name' ,
    'n',
    1,
    'character',
    'outpath'      ,
    'o',
    1,
    'character',
    'k'            ,
    'k',
    1,
    "integer",
    'min_size'     ,
    'm',
    1,
    "integer",
    'n_boot'       ,
    'b',
    1,
    "integer"
  ),
  byrow = TRUE,
  ncol = 4
)

opt = getopt(spec)



meta_cell_func <-
  function(single_cell_fn,
           genes_fn,
           K,
           K2,
           min_mc_size,
           n_boot,
           dataset_name,
           outpath) {
    message(dataset_name)
    message(K)
    message(K2)
    message(min_mc_size)
    message(n_boot)
    
    #Create loom file connection
    loom.data <- connect(single_cell_fn, mode = 'r')
    message("Opened Loom Connection")
    seurat.obj.from.loom <-
      as.Seurat(loom.data, cells = 'obs_names', features = 'var_names')
    message("Created Seurat Object")
    
    genes <- read.table(genes_fn, stringsAsFactors = F)[, 1]
    
    sce <- as.SingleCellExperiment(seurat.obj.from.loom)
    message("Created Single Cell Experiment")
    #Subset to gene names that I want to use and Neurons only
    sce <-
      sce[genes, sce@colData$class_label %in% c("GABAergic", "Glutamatergic")]
    mat <- scm_import_sce_to_mat(sce)
    message("Created tgSCmat")
    
    #Initialize Outpath for data and figures
    scdb_init(outpath, force_reinit = T)
    scfigs_init(outpath)
    
    #Put expression data in outpath
    scdb_add_mat(dataset_name, mat)
    
    
    mcell_plot_umis_per_cell(dataset_name)
    message("Plotted UMIs per cell")
    
    #Prepare GeneSet (all genes)
    mcell_add_gene_stat(gstat_id = "all_genes",
                        mat_id = dataset_name,
                        force = T)
    mcell_gset_filter_varmean(
      gstat_id = "all_genes",
      gset_id = "all_geneset",
      T_vm = 0,
      force_new = T
    )
    mcell_plot_gstats(gstat_id = "all_genes", gset_id = "all_geneset")
    message("Plotting gene statistics")
    
    #Creates Initial Graph
    mcell_add_cgraph_from_mat_bknn(
      mat_id = dataset_name,
      gset_id = "all_geneset",
      graph_id = paste(dataset_name, "graph", sep =
                         '_'),
      K = K
    )
    message("Created Initial Graph")
    #Resample Graph on subset to get concensus
    mcell_coclust_from_graph_resamp(
      coc_id = paste(dataset_name, "coclust", sep = '_'),
      graph_id = paste(dataset_name, "graph", sep =
                         '_'),
      min_mc_size = min_mc_size,
      n_resamp = n_boot,
      p_resamp = .75
    )
    message("Completed Bootstrapping")
    mcell_mc_from_coclust_balanced(
      coc_id = paste(dataset_name, "coclust", sep = '_'),
      mat_id = dataset_name,
      mc_id = paste(dataset_name, "metacells", sep =
                      '_'),
      K = K2,
      min_mc_size = min_mc_size
    )
    
    message("Done")
    message(dataset_name)
    load(paste(outpath, 'mc.', dataset_name, '_metacells.Rda', sep = ''))
    object@mc_fp %>% write.csv(
      paste(outpath, dataset_name, 'meta_cell_matrix.csv', sep =
              ''),
      row.names = T,
      quote = F
    )
    object@mc %>% write.csv(
      paste(outpath, dataset_name, 'meta_cell_assignment.csv', sep =
              ''),
      row.names = T,
      quote = F
    )
  }
meta_cell_func(opt$loomfile, opt$genefile, opt$k, opt$k, opt$min_size, opt$n_boot, opt$dataset_name, opt$outpath)