
source("datasets.R")
source("differential_expression.R")


main = function() {
    dir.create("M1", showWarnings = FALSE)
    make_meta_lists(biccn_datasets(), "M1")
}
                                 
make_meta_lists = function(dataset_names, output_dir) {
    clusters = make.names(c("GABAergic", "Glutamatergic", "Non-Neuronal",
                            excitatory_subclasses(), inhibitory_subclasses()))
    for (c in clusters) {
        file_list = paste0(dataset_names, "/", c, ".txt")
        file_list = file_list[file.exists(file_list)]
        create_export_list_de(file_list, c, output_dir)
        create_export_list_fc(file_list, c, output_dir)
    }
}
                                 
create_export_list_de = function(file_list, file_suffix, meta_dir) {
    result = create_meta_gene_list_de(file_list)
    export_meta_list(result, file.path(meta_dir, paste0(file_suffix, "_markers_de.txt")))
}
                                 
create_export_list_fc = function(file_list, file_suffix, meta_dir) {
    result = create_meta_gene_list_fc(file_list)
    export_meta_list(result, file.path(meta_dir, paste0(file_suffix, "_markers_fc.txt")))
}
                                 
if (!interactive()) {
    main()
}