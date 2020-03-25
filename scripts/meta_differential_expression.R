


aggregate_de_stats = function(file_list) {
    stats = lapply(file_list, read.table, header = TRUE)
    names(stats) = sapply(file_list, basename)
    
    aurocs = extract_stat(stats, "auroc")
    gene_occurrence = compute_feature_occurrence(aurocs)
    auroc_matrix = list_to_matrix(aurocs)
    
    return(list(
        auroc_matrix = auroc_matrix,
        fc_matrix = list_to_matrix(extract_stat(stats, "fold_change")),
        fdr_matrix = list_to_matrix(extract_stat(stats, "log_fdr")),
        detection_matrix = list_to_matrix(extract_stat(stats, "detection_rate")),
        expression_matrix = list_to_matrix(extract_stat(stats, "average_expression")),
        precision_matrix = list_to_matrix(extract_stat(stats, "precision")),
        gene_occurrence = gene_occurrence[rownames(auroc_matrix)]
    ))
}

extract_stat = function(stats, stat_name) {
    result = lapply(stats, function(s) {
        res = s[, stat_name]
        names(res) = as.character(s$gene)
        return(res)
    })
    return(result)
}

compute_feature_occurrence = function(stat_list) {
    return(c(table(unlist(lapply(stat_list, names)))))
}

list_to_matrix = function(stat_list) {
    features = Reduce(union, lapply(stat_list, names))
    result = matrix(0, nrow = length(features), ncol = length(stat_list),
                    dimnames = list(features, names(stat_list)))
    for (i in seq_along(stat_list)) {
        result[names(stat_list[[i]]), i] = stat_list[[i]]
    }
    return(result)
}

summarize_de_stats = function(aggregate_stats, recurrence_matrix, order_by = "average_auroc") {
    result = data.frame(
        gene = rownames(aggregate_stats$auroc_matrix),
        recurrence = rowSums(recurrence_matrix),
        average_auroc = rowSums(aggregate_stats$auroc_matrix) / aggregate_stats$gene_occurrence,
        average_fold_change = rowSums(aggregate_stats$fc_matrix) / aggregate_stats$gene_occurrence,
        average_detection = rowSums(aggregate_stats$detection_matrix) / aggregate_stats$gene_occurrence,
        average_expression = rowSums(aggregate_stats$expression_matrix) / aggregate_stats$gene_occurrence,
        average_precision = rowSums(aggregate_stats$precision_matrix) / aggregate_stats$gene_occurrence,
        maximal_recurrence = aggregate_stats$gene_occurrence
    )
    result = cbind(result, recurrence_matrix)
    result = result[order(result$recurrence, result[[order_by]], decreasing = TRUE),]
    return(result)
}

create_meta_gene_list_de = function(file_list, fc_threshold = 4, fdr_threshold = 0.05) {
    aggregate_stats = aggregate_de_stats(file_list)
    recurrence_matrix = compute_de_recurrence(aggregate_stats, fc_threshold, fdr_threshold)
    return(list(meta_list = summarize_de_stats(aggregate_stats, recurrence_matrix),
                files = sapply(file_list, basename),
                method = paste0("recurrence among DE genes (FC > ", fc_threshold, ", FDR < ", fdr_threshold,")")))
}

compute_de_recurrence = function(aggregate_stats, fc_threshold = 4, fdr_threshold = 0.05) {
    result = 1*(aggregate_stats$fc_matrix > fc_threshold & aggregate_stats$fdr_matrix < log(fdr_threshold))
    dimnames(result) = dimnames(aggregate_stats$fc_matrix)
    return(result)
}

create_meta_gene_list_fc = function(file_list, fc_threshold = 4, detection_threshold = 0.05) {
    aggregate_stats = aggregate_de_stats(file_list)
    recurrence_matrix = compute_fc_recurrence(aggregate_stats, fc_threshold, detection_threshold)
    return(list(meta_list = summarize_de_stats(aggregate_stats, recurrence_matrix, "average_fold_change"),
                files = sapply(file_list, basename),
                method = paste0("recurrence among DE genes (FC > ", fc_threshold, ", detection_rate > ", detection_threshold, ")")))
}

compute_fc_recurrence = function(aggregate_stats, fc_threshold = 4, detection_threshold = 0.05) {
    result = 1*(aggregate_stats$fc_matrix > fc_threshold & aggregate_stats$detection_matrix > detection_threshold)
    dimnames(result) = dimnames(aggregate_stats$fc_matrix)
    return(result)
}

export_meta_list = function(meta_list, filename) {
    write_meta_header(meta_list, filename)
    write.table(meta_list$meta_list, filename, append = TRUE, row.names = FALSE)
}

write_meta_header = function(meta_list, filename) {
    lines = paste0(
        "# Marker list generated on ", Sys.Date(), ". ",
        "Ordered by ", meta_list$method, ". ",
        "Based on: ", paste(meta_list$files, collapse = ", "), "."
    )
    writeLines(lines, filename)
}
