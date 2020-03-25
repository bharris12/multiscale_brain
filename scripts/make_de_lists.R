
library(tidyverse)

source("datasets.R")
source("differential_expression.R")


main = function() {
    make_de_lists_biccn()
}

make_de_lists_biccn = function() {
    dataset_names = set_names(biccn_datasets())
    for (name in dataset_names) {
        print(name)
        dir.create(name, showWarnings = FALSE)
        gc()
        dataset = load_biccn_dataset(name)
        gc()
        make_de_lists_dataset(dataset, name)
    }
}

make_de_lists_dataset = function(dataset, name) {
    make_de_list_class(dataset, name)
    gc()
    make_de_list_subclass(dataset, name, "Glutamatergic")
    gc()
    make_de_list_subclass(dataset, name, "GABAergic")
}

make_de_list_class = function(dataset, name) {
    label_matrix = design_matrix(as.character(dataset$class_label))
    de_stats = compute_de_blockwise(dataset, label_matrix)
    export_de_stats(de_stats, make.names(colnames(label_matrix)), paste0(name, "/"))
}

make_de_list_subclass = function(dataset, name, class_name) {
    dataset = dataset[, dataset$class_label == class_name]
    label_matrix = design_matrix(dataset$common_subclass)
    de_stats = compute_de_blockwise(dataset, label_matrix)
    export_de_stats(de_stats, make.names(colnames(label_matrix)), paste0(name, "/"))
}
                                 
if (!interactive()) {
    main()
}