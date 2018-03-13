##' Measuring mutually exclusive transcriptional phenotypes
##' @rdname measureTrPhe
##' @name measureTrPhe
##' @author Vitalii Kleshchevnikov
##' @description \code{measureTrPhe} idenfifies mutually exclusive transcriptional phenotypes in single cell data (normalised read counts, SingleCellExperiment)
##' @param data object of SingleCellExperiment containing single cell data (normalised read counts, cells already assigned to clusters)
##' @param method method for detecting differentially expressed genes. Currently only Wilcox test is implemented \link[stats]{wilcox.test}.
##' @param mode compare clusters to each other (\code{"pairwise"}) or \code{"one_vs_all"}
##' @param cutoff FDR-corrected p-value cutoff
##' @param assays_matrix_name name of the matrix in \code{assays(data)} that stores normalised read counts
##' @param colData_cluster_col name of the column in \code{colData(data)} that stores cluster assignment of cells
##' @param pval_corr_method multiple hypothesis p-value correction method. Details: \link[stats]{p.adjust} - method.
##' @param low_exprs_threshold threshold (normalised read count) below which gene doesn't qualify as being expressed in a cell
##' @param low_exprs_cells remove genes that are expressed at a level higher or equal \code{low_exprs_threshold} in fewer than \code{low_exprs_cells} cells at each between-cluster comparison
##' @return data.table containing the perturbation details from the Connectivity map project. Only one overexpression experiment per gene and condition is retained.
##' @import data.table
##' @import SingleCellExperiment
##' @import parallel
##' @export measureTrPhe
##' @examples
##' library(ArrayExpress)
##' library(SingleCellExperiment)
##' library(data.table)
##' file_paths = getAE("E-MTAB-6153", type = "processed", path = "../regulatory_networks_by_cmap/data/organogenesis_scRNAseq", local = T)
##' # keep only normalised counts
##' file_paths$processedFiles = file_paths$processedFiles[file_paths$processedFiles == "normalisedCounts.tsv"]
##' # get the list of column names
##' cnames = getcolproc(file_paths)
##' data = readEMTAB6153ProcData(path = "../regulatory_networks_by_cmap/data/organogenesis_scRNAseq", procFile = "normalisedCounts.tsv", procol = cnames)
measureTrPhe = function(data, method = "wilcox", mode = c("pairwise", "one_vs_all")[1], cutoff = 0.05, assays_matrix_name = "norm_counts", colData_cluster_col = "clusters", pval_corr_method = "fdr", low_exprs_threshold = 0.1, low_exprs_cells = 6, cutoff2nd = 0.5){
  clusters = colData(data)[,colData_cluster_col]
  combinations = clusterCOMBS(clusters, mode)
  norm_counts = assays(data)[[assays_matrix_name]]


  measureTrPheSingle(norm_counts, combinations, clusters, combinations_ind = 1,
                     cutoff = cutoff, pval_corr_method = pval_corr_method,
                     low_exprs_threshold = low_exprs_threshold, low_exprs_cells = low_exprs_cells,
                     method = method, cutoff2nd = cutoff2nd)


  microbenchmark::microbenchmark({p1 <- measureTrPheSingle(norm_counts, combinations, clusters, combinations_ind = 1,
                                                           cutoff = cutoff, pval_corr_method = pval_corr_method,
                                                           low_exprs_threshold = low_exprs_threshold, low_exprs_cells = low_exprs_cells,
                                                           method = method, cutoff2nd = cutoff2nd)},
                                 {p2 <- measureTrPheSingle(norm_counts, combinations, clusters, combinations_ind = 10,
                                                           cutoff = cutoff, pval_corr_method = pval_corr_method,
                                                           low_exprs_threshold = low_exprs_threshold, low_exprs_cells = low_exprs_cells,
                                                           method = method, cutoff2nd = cutoff2nd)},
                                 {p3 <- measureTrPheSingle(norm_counts, combinations, clusters, combinations_ind = 100,
                                                           cutoff = cutoff, pval_corr_method = pval_corr_method,
                                                           low_exprs_threshold = low_exprs_threshold, low_exprs_cells = low_exprs_cells,
                                                           method = method, cutoff2nd = cutoff2nd)},
                                 {p4 <- measureTrPheSingle(norm_counts, combinations, clusters, combinations_ind = 200,
                                                           cutoff = cutoff, pval_corr_method = pval_corr_method,
                                                           low_exprs_threshold = low_exprs_threshold, low_exprs_cells = low_exprs_cells,
                                                           method = method, cutoff2nd = cutoff2nd)}, times = 5)

}

##' @rdname measureTrPhe
##' @name measureTrPheSingle
##' @description \code{measureTrPheSingle}: idenfifies mutually exclusive transcriptional phenotypes for 2 clusters of single cells
##' @param clusters character vector, clusters (cell types) present in the data
##' @param mode
##' @return \code{measureTrPheSingle}: data.table containing phenotype id (phenotypes) and which genes are assigned to them
##' @import data.table
##' @export measureTrPheSingle
measureTrPheSingle = function(norm_counts, combinations, clusters, combinations_ind = 1, cutoff = 0.05, pval_corr_method = "fdr", low_exprs_threshold = 1, low_exprs_cells = 6, method = c("wilcox","ks")[1], cutoff2nd = 0.5){
  if(length(clusters) != ncol(norm_counts)) stop("measureTrPheSingle: length of vector \"clusters\" should be equal to the number of columns in the read count matrix")
  #norm_counts_subset = norm_counts[, clusters %in% c(combinations[combinations_ind, phenotypes_cluster1],
  #                                                   combinations[combinations_ind, phenotypes_cluster2])]

  norm_counts_subset = norm_counts[, grepl(combinations[combinations_ind, phenotypes_cluster1],
                                           colnames(norm_counts)) |
                                     grepl(combinations[combinations_ind, phenotypes_cluster2],
                                           colnames(norm_counts))]
  # filter out genes that are not expressed in both of these clusters or are expressed in a few cells
  rowsums = rowSums2(norm_counts_subset >= low_exprs_threshold)
  norm_counts_subset = norm_counts_subset[rowsums > low_exprs_cells,]
  if(method == "wilcox"){
    #pVals <- apply(norm_counts_subset,
    #               2, function(x){
    #                 wilcox.test( x[clusters == combinations[combinations_ind, phenotypes_cluster1]],
    #                              x[clusters == combinations[combinations_ind, phenotypes_cluster2]],
    #                              alternative = "greater")$p.value} )
    pVals <- apply(norm_counts_subset,
                   1, function(x){
                     wilcox.test( x[grepl(combinations[combinations_ind, phenotypes_cluster1],
                                          colnames(norm_counts_subset))],
                                  x[grepl(combinations[combinations_ind, phenotypes_cluster2],
                                          colnames(norm_counts_subset))],
                                  alternative = "greater")$p.value} )
  }
  if(method == "ks"){
    pVals <- apply(norm_counts_subset,
                   1, function(x){
                     ks.test( x[grepl(combinations[combinations_ind, phenotypes_cluster1],
                                      colnames(norm_counts_subset))],
                              x[grepl(combinations[combinations_ind, phenotypes_cluster2],
                                      colnames(norm_counts_subset))],
                              alternative = "less")$p.value} )
  }

  # multiple testing correction - do for all 380 comparisons
  pVals_corr <- p.adjust(pVals, method = pval_corr_method)
  passed = pVals_corr < cutoff
  if(sum(passed) == 0){
    warning(paste0(combinations[combinations_ind, phenotypes], " has no genes below fdr threshold of ", cutoff))
    #warning(paste0(combinations[combinations_ind, phenotypes], " has no genes below fdr threshold of ", cutoff,", choosing genes using uncorrected p-value < ",cutoff2nd))
    passed = pVals < cutoff2nd
    if(sum(passed) > 0){
      phenotype2gene = data.table(phenotypes = combinations[combinations_ind, phenotypes],
                                  genes = NA,# names(pVals)[passed],
                                  pVals = NA)# pVals[passed])
    } else {
      phenotype2gene = data.table(phenotypes = combinations[combinations_ind, phenotypes],
                                  genes = NA,
                                  pVals = NA)
    }
  } else {
    phenotype2gene = data.table(phenotypes = combinations[combinations_ind, phenotypes],
                                genes = names(pVals_corr)[passed],
                                pVals = pVals_corr[passed])
  }
  phenotype2gene
}

##' @rdname measureTrPhe
##' @name clusterCOMBS
##' @description \code{clusterCOMBS}: produces a data.table specifying all combinations of clusters for differential expression analysis.
##' @param clusters character vector, clusters (cell types) present in the data
##' @param mode
##' @return \code{clusterCOMBS}: data.table containing phenotype id (phenotypes) and which clusters are to be compared
##' @import data.table
##' @export clusterCOMBS
clusterCOMBS = function(clusters, mode = c("pairwise", "one_vs_all")[1]) {
  clusters = unique(clusters)
  if(mode == "one_vs_all") clusters = c(clusters, "all")
  combinations = character()
  for (cluster1 in clusters) {
    for (cluster2 in clusters) {
      if(cluster1 != cluster2) {
        combinations = c(combinations, paste0(cluster1, "_NOT", cluster2))
      }
    }
  }
  if(mode == "one_vs_all") combinations = grep("all", combinations, value = T)
  combinations = data.table(phenotypes = combinations)
  #combinations[, iterations_cluster1 := {
  #  temp = unlist(tstrsplit(phenotypes, "_NOT"))
  #  temp = temp[order(temp)][1]
  #}, by = phenotypes]
  #combinations[, iterations_cluster2 := {
  #  temp = unlist(tstrsplit(phenotypes, "_NOT"))
  #  temp = temp[order(temp)][2]
  #}, by = phenotypes]
  combinations[, phenotypes_cluster1 := unlist(tstrsplit(phenotypes, "_NOT"))[1], by = phenotypes]
  combinations[, phenotypes_cluster2 := unlist(tstrsplit(phenotypes, "_NOT"))[2], by = phenotypes]
  # setorder(combinations, iterations_cluster1, iterations_cluster2)
  combinations
}
