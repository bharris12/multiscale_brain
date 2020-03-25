
source("aurocs.R")
source("meta_differential_expression.R")


compute_de = function(expression, design_matrix, two_tailed = TRUE, tie_correction = FALSE) {
    aurocs = compute_aurocs(t(expression), design_matrix, return_tie_correction = tie_correction)
    if (tie_correction) {
        p_values = t(auroc_p_value(aurocs$aurocs, design_matrix, two_tailed, tie_correction = aurocs$tie_correction, log.p = TRUE))
        aurocs = aurocs$aurocs
    } else {
        p_values = t(auroc_p_value(aurocs, design_matrix, two_tailed, log.p = TRUE))
    }
    
    population_size = colSums(design_matrix)

    average = average_expression(expression, design_matrix)
    uncentered_var = average_expression(expression**2, design_matrix)
    standard_error = sqrt((uncentered_var$positives - average$positives**2) / population_size)
    
    binary_expression = expression > 0
    m_binary = average_expression(binary_expression, design_matrix)
    n_expressing_cells = binary_expression %*% design_matrix
    binary_precision = n_expressing_cells / rowSums(binary_expression)
    binary_recall = t(t(n_expressing_cells) / population_size)
      
    return(list(
        gene = rownames(expression),
        population_size = population_size,
        population_fraction = population_size / sum(population_size),
        auroc = t(aurocs),
        log_p_value = p_values,
        log_fdr = matrix(my_fdr(p_values, log.p = TRUE), nrow = nrow(p_values), dimnames = dimnames(p_values)),
        average_expression = average$positives,
        se_expression = standard_error,
        fold_change = (average$positives+1) / (average$negatives+1),
        detection_rate = m_binary$positives,
        fold_change_detection = (m_binary$positives + 0.001) / (m_binary$negatives + 0.001),
        precision = binary_precision,
        recall = binary_recall
    ))
}

average_expression = function(expression, design_matrix) {
    positive_expression = expression %*% design_matrix
    negative_expression = rowSums(expression) - positive_expression
    return(list(
        positives = scale(positive_expression, center=FALSE, scale=colSums(design_matrix)),
        negatives = scale(negative_expression, center=FALSE, scale=colSums(1-design_matrix))
    ))
}

# adpated from p.adjust
my_fdr = function(p_values, log.p = FALSE) {
    i = length(p_values):1L
    o = order(p_values, decreasing = TRUE)
    n = length(p_values)
    result = rep(0, n)
    if (log.p) {
        result[o] = pmin(0, cummin(log(n) - log(i) + p_values[o]))
    } else {
        result[o] = pmin(1, cummin(n/i * p_values[o]))
    }
    return(result)
}

compute_de_blockwise = function(dataset, labels, two_tailed = TRUE, tie_correction = FALSE) {
    blocks = unique(c(seq(1, nrow(dataset), floor(2e9 / ncol(dataset))), nrow(dataset)))
    result = lapply(seq_len(length(blocks) - 1), function(i) compute_de(
        as.matrix(assay(dataset[blocks[i]:blocks[i+1],], "cpm")), labels, two_tailed, tie_correction)
    )
    return(fuse_de_stats(result))
}

fuse_de_stats = function(stat_list) {
    result = lapply(seq_along(stat_list[[1]]), function(i) {
        Reduce(rbind, lapply(stat_list, function(s) as.matrix(s[[i]])))
    })
    names(result) = names(stat_list[[1]])
    return(result)
}
               
export_de_stats = function(de_stats, population_names, file_prefix) {
    sapply(seq_along(population_names), function(i) export_population(de_stats, i, population_names, file_prefix))
    return(NULL)
}

export_population = function(de_stats, index, population_names, file_prefix) {
    filename = paste0(file_prefix, population_names[index], ".txt")
    data_ = get_population(de_stats, index)
    write_header(population_names[index], population_names[-index], filename)
    write.table(data_, filename, append = TRUE, row.names = FALSE)
}
                              
get_population = function(de_stats, index) {
    population_stats = lapply(4:length(de_stats), function(i) de_stats[[i]][, index])
    names(population_stats) = names(de_stats)[-(1:3)]
    populations_stats = as.data.frame(population_stats)
    data_ = cbind(data.frame(gene = de_stats$gene,
                             population_size = de_stats$population_size[index],
                             population_fraction = de_stats$population_fraction[index]),
                  population_stats)
    data_ = data_[order(data_$auroc, decreasing = TRUE), ]
    return(data_)
}
           
write_header = function(reference, outgroups, filename) {
    lines = paste0(
        "# Differential expression statistics generated on ", Sys.Date(), ". ",
        "Reference: ", reference, ", outgroups: ", paste(outgroups, collapse = ", "), "."
    )
    writeLines(lines, filename)
}
