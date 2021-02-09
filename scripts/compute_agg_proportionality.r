library(tidyverse)
library(SingleCellExperiment)
library(propr)
library(argparser)
library(testit)

source('~/r_utils/read_loom.R')
source('~/r_utils/rank.r')


options(echo = T)
parser <-
    arg_parser(description = 'Dataset Name and Aggregation Level')
parser <- add_argument(parser,
                       '--dataset',
                       help = 'Name of dataset')
parser <-
    add_argument(parser,
                 '--level',
                 help = 'Level of aggregation, either class, subclass, or joint_cluster')
args <- parse_args(parser)
dataset <-
    as.vector(args$dataset) %>% substr(1, nchar(args$dataset) - 12)
dataset <- strsplit(dataset, '/')[[1]][7]
agg_level <- args$level

##Constants

genes_fn <-
    '/data/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv'


##Load Data
genes <- read.csv(genes_fn, stringsAsFactors = F, header=F)[, 1]

dataset_fn <- args$dataset
sce <- get_sce_alt(dataset_fn)
sce <- sce[, sce$class_label %in% c('GABAergic', 'Glutamatergic')]
sce <- sce[genes, ]

agg_levels <- unique(colData(sce)[, agg_level])
agg_nw <-
    matrix(
        0,
        nrow = length(genes),
        ncol = length(genes),
        dimnames = list(genes, genes)
    )
for (subset in agg_levels) {
    prop <-
        propr(as.matrix(t(counts(sce[, colData(sce)[, agg_level] == subset]))), p =
                  0)
    nw <- rank_matrix(prop@matrix)
    agg_nw <- agg_nw + nw
    
    rm(prop)
    rm(nw)
    gc()
}
agg_nw <- rank_matrix(agg_nw)
output_fn <-
    paste0(
        '~/biccn_paper/data/networks/proportionality/proportionality_nw_',
        dataset,
        '_',
        agg_level,
        '.csv.gz'
    )
write.csv(agg_nw, gzfile(output_fn))