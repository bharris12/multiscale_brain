library(gemmaAPI)
library(tidyverse)

datasets<-read.table('~/biccn_paper/data/bulk_rna/selected_gemma_datasets.txt',stringsAsFactors = F)[,1]
outpath<-'~/biccn_paper/data/bulk_rna/gemma_data/'
message(outpath)

datasets %>% sapply(function(x){datasetInfo(x,request= 'data',return= FALSE,memoised=TRUE, file = paste0(outpath,x))})

