
library(SingleCellExperiment)


biccn_datasets = function() {
    c("zeng_smart_cells", "zeng_smart_nuclei", "zeng_10x_cells_v2", "zeng_10x_nuclei_v2",
      "macosko_10x_nuclei_v3", "zeng_10x_cells_v3", "zeng_10x_nuclei_v3")
}

load_biccn_dataset = function(dataset_name) {
    dataset = readRDS(paste0("~/data/biccn/parsed_data/", dataset_name, ".rds"))
    dataset = dataset[, dataset$class_label %in% c("GABAergic", "Glutamatergic", "Non-Neuronal")]
    assay(dataset, "cpm") = convert_to_cpm(assay(dataset))
    labels = dataset$subclass_label
    labels[dataset$class_label == "Non-Neuronal"] = "Non-Neuronal"
    labels[labels == "L5 NP"] = "L5/6 NP"
    dataset$common_subclass = labels
    dataset$class = dataset$class_label
    dataset$subclass = dataset$subclass_label
    dataset$cluster = dataset$cluster_label
    return(dataset)
}

convert_to_cpm = function(M) {
    normalization_factor = Matrix::colSums(M) / 1000000
    if (class(M) == "dgCMatrix") {
        M@x = M@x / rep.int(normalization_factor, diff(M@p))
        return(M)
    } else {
        return(scale(M, center = FALSE, scale = normalization_factor))
    }
}

excitatory_subclasses = function() {
    return(c("L2/3 IT", "L5 IT", "L5/6 NP", "L5 ET", "L6 CT", "L6 IT", "L6 IT Car3", "L6b"))
}

inhibitory_subclasses = function() {
    return(c("Lamp5", "Pvalb", "Sst", "Vip", "Sncg"))
}

non_neuronal_subclasses = function() {
    return(c("Astro", "Oligo", "OPC", "Micro", "Endo", "SMC", "VLMC", "Macrophage", "Peri"))
}

