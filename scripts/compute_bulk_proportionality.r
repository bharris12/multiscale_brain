library(tidyverse)
library(SingleCellExperiment)
library(propr)
library(argparser)
library(testit)


source('~/r_utils/rank.r')

genes_fn <-
    '/data/bharris/biccn_paper/data/highly_expressed_7_datasets_75k.csv'
genes <- read.csv(genes_fn, stringsAsFactors = F, header = F)[, 1]

datasets_used <-
    read.csv(
        '~/biccn_paper/data/bulk_rna/datasets_used.csv',
        stringsAsFactors = F,
        header = F
    )[, 1]

agg_nw <-
    matrix(
        0,
        nrow = length(genes),
        ncol = length(genes),
        dimnames = list(genes, genes)
    )

compute_bulk_prop <- function(dataset) {
    expr <- read.csv(
        paste0(
            '~/biccn_paper/data/bulk_rna/gemma_genes_reduced/',
            dataset,
            '_expression.csv.gz'
        ),
        row.names = 1
    )
    expr <- expr + abs(min(expr))
    prop <- propr(t(expr), p = 0)@matrix
    rank_matrix(prop)
}

for (dataset in datasets_used) {
    agg_nw <- agg_nw + compute_bulk_prop(dataset)
}
agg_nw <- rank_matrix(agg_nw)
write.csv(
    agg_nw,
    gzfile(
        '~/biccn_paper/data/networks/proportionality/bulk_agg_nw.csv.gz'
    )
)