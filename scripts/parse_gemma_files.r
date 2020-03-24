library(gemmaAPI)
library(tidyverse)

datasets <-
  read.table('~/biccn_paper/data/bulk_rna/selected_gemma_datasets.txt',
             stringsAsFactors = F)[, 1]

for (dataset in datasets) {
  message(dataset)
  full_dataset <-
    readExpression(paste0('~/biccn_paper/data/bulk_rna/gemma_data/', dataset))
  genes <- full_dataset$GeneSymbol
  samples_slice <- 7:dim(full_dataset)[2]
  expression <- as.matrix(full_dataset[, samples_slice])
  rownames(expression) <- genes
  write.csv(
    expression,
    paste0(
      '~/biccn_paper/data/bulk_rna/gemma_parsed_expression/',
      dataset,
      '_expression.csv'
    ),
    quote = F
  )
}
