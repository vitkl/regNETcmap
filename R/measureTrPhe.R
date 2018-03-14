##' Measuring mutually exclusive transcriptional phenotypes
##' @rdname measureTrPhe
##' @name measureTrPhe
##' @author Vitalii Kleshchevnikov
##' @description \code{measureTrPhe} idenfifies mutually exclusive transcriptional phenotypes in single cell data (normalised read counts, SingleCellExperiment): genes expressed in one cell type but not the other (or all other) using one-tailed Wilcox, KS or other differential expression tests
##' @param data object of SingleCellExperiment containing single cell data (normalised read counts, cells already assigned to clusters)
##' @param method method for detecting differentially expressed genes. Currently only Wilcox and KS tests are implemented \link[stats]{wilcox.test}, \link[stats]{ks.test}.
##' @param mode compare clusters to each other (\code{"pairwise"}) or \code{"one_vs_all"}
##' @param cutoff FDR-corrected p-value cutoff
##' @param assays_matrix_name name of the matrix in \code{assays(data)} that stores normalised read counts
##' @param colData_cluster_col name of the column in \code{colData(data)} that stores cluster assignment of cells
##' @param pval_corr_method multiple hypothesis p-value correction method. Details: \link[stats]{p.adjust} - method.
##' @param low_exprs_threshold threshold (normalised read count) below which gene doesn't qualify as being expressed in a cell
##' @param low_exprs_cells remove genes that are expressed at a level higher or equal \code{low_exprs_threshold} in fewer than \code{low_exprs_cells} cells at each between-cluster comparison
##' @param n_cores number of cores to be used in parallel processing (over combinations of clusters). More details: \link[parallel]{parLapply}, \link[parallel]{makeCluster}, \link[parallel]{detectCores}
##' @return data.table containing the perturbation details from the Connectivity map project. Only one overexpression experiment per gene and condition is retained.
##' @import data.table
##' @import SingleCellExperiment
##' @import parallel
##' @export measureTrPhe
##' @examples
##' library(ArrayExpress)
##' library(SingleCellExperiment)
##' library(data.table)
##' library(ggplot2)
##' file_paths = getAE("E-MTAB-6153", type = "processed", path = "../regulatory_networks_by_cmap/data/organogenesis_scRNAseq", local = T)
##' # keep only normalised counts
##' file_paths$processedFiles = file_paths$processedFiles[file_paths$processedFiles == "normalisedCounts.tsv"]
##' # get the list of column names
##' cnames = getcolproc(file_paths)
##' data = readEMTAB6153ProcData(path = "../regulatory_networks_by_cmap/data/organogenesis_scRNAseq", procFile = "normalisedCounts.tsv", procol = cnames)
##' pVals = measureTrPhe(data, method = "wilcox", mode = c("pairwise", "one_vs_all")[1], cutoff = 1, assays_matrix_name = "norm_counts", colData_cluster_col = "clusters", pval_corr_method = "fdr", low_exprs_threshold = 0.1, low_exprs_cells = 6, n_cores = detectCores() - 1)
##' qplot(x = pVals$diff_median, y = -log10(pVals$pVals), geom = "bin2d", xlim = c(-1,50), ylim = c(-1, 300), bins = 150) + theme_light()
measureTrPhe = function(data, method = "wilcox", mode = c("pairwise", "one_vs_all")[1], cutoff = 0.05, assays_matrix_name = "norm_counts", colData_cluster_col = "clusters", pval_corr_method = "fdr", low_exprs_threshold = 0.1, low_exprs_cells = 6, n_cores = detectCores() - 1){
  clusters = colData(data)[,colData_cluster_col]
  combinations = clusterCOMBS(clusters, mode)
  norm_counts = assays(data)[[assays_matrix_name]]

  # set up parallel processing
  # create cluster
  cl <- makeCluster(n_cores)
  # get library support needed to run the code
  clusterEvalQ(cl, {library(regNETcmap)})
  # put objects in place that might be needed for the code
  clusterExport(cl, c("norm_counts", "combinations", "cutoff", "pval_corr_method", "low_exprs_threshold", "low_exprs_cells", "method"), envir=environment())
  pVals = parLapply(cl, X = 1:nrow(combinations),
                    fun = function(ind) {
                      measureTrPheSingle(norm_counts = norm_counts, combinations = combinations,
                                         combinations_ind = ind,
                                         cutoff = cutoff, pval_corr_method = pval_corr_method,
                                         low_exprs_threshold = low_exprs_threshold,
                                         low_exprs_cells = low_exprs_cells,
                                         method = method)
                    })
  # stop cluster
  stopCluster(cl)

  # combine results
  pVals = Reduce(rbind, pVals)
  # multiple testing correction - do for all 380 comparisons
  pVals[, pVals_corr := p.adjust(pVals, method = pval_corr_method)]
  pVals[, passed := pVals_corr < cutoff]
  pVals
}

##' @rdname measureTrPhe
##' @name measureTrPheSingle
##' @description \code{measureTrPheSingle}: idenfifies mutually exclusive transcriptional phenotypes for 2 clusters of single cells
##' @param clusters character vector, clusters (cell types) present in the data
##' @param combinations data.table containing phenotype id (phenotypes) and which clusters are to be compared
##' @param combinations_ind which combination should be analysed?
##' @return \code{measureTrPheSingle}: data.table containing phenotype id (phenotypes) and which genes are assigned to them
##' @import data.table
##' @export measureTrPheSingle
measureTrPheSingle = function(norm_counts, combinations, combinations_ind = 1, cutoff = 0.05, pval_corr_method = "fdr", low_exprs_threshold = 1, low_exprs_cells = 6, method = c("wilcox","ks")[1]){

  if(combinations[combinations_ind, phenotypes_cluster1] == "all"){
    phenotypes_cluster2 = combinations[combinations_ind, phenotypes_cluster2]
    norm_counts_subset = norm_counts
  } else if(combinations[combinations_ind, phenotypes_cluster2] == "all"){
    phenotypes_cluster1 = combinations[combinations_ind, phenotypes_cluster1]
    norm_counts_subset = norm_counts
  } else {
    phenotypes_cluster1 = combinations[combinations_ind, phenotypes_cluster1]
    phenotypes_cluster2 = combinations[combinations_ind, phenotypes_cluster2]
    norm_counts_subset = norm_counts[, grepl(phenotypes_cluster1, colnames(norm_counts)) |
                                       grepl(phenotypes_cluster2, colnames(norm_counts))]
  }

  # filter out genes that are not expressed in both of these clusters or are expressed in a few cells
  rowsums = matrixStats::rowSums2(norm_counts_subset >= low_exprs_threshold)
  norm_counts_subset = norm_counts_subset[rowsums > low_exprs_cells,]
  if(method == "wilcox"){
    if(combinations[combinations_ind, phenotypes_cluster1] == "all"){
      pVals <- apply(norm_counts_subset,
                     1, function(mat){
                       x = mat[!grepl(phenotypes_cluster2, colnames(norm_counts_subset))]
                       y = mat[grepl(phenotypes_cluster2, colnames(norm_counts_subset))]
                       w = wilcox.test(x, y, alternative = "greater")
                       c(w$p.value, w$statistic, median(x) - median(y))} )
    } else if(combinations[combinations_ind, phenotypes_cluster2] == "all"){
      pVals <- apply(norm_counts_subset,
                     1, function(mat){
                       x = mat[grepl(phenotypes_cluster1, colnames(norm_counts_subset))]
                       y = mat[!grepl(phenotypes_cluster1, colnames(norm_counts_subset))]
                       w = wilcox.test(x, y, alternative = "greater")
                       c(w$p.value, w$statistic, median(x) - median(y))} )
    } else {
      pVals <- apply(norm_counts_subset,
                     1, function(mat){
                       x = mat[grepl(phenotypes_cluster1, colnames(norm_counts_subset))]
                       y = mat[grepl(phenotypes_cluster2, colnames(norm_counts_subset))]
                       w = wilcox.test(x, y, alternative = "greater")
                       c(w$p.value, w$statistic, median(x) - median(y))} )
    }
  }
  if(method == "ks"){
    if(combinations[combinations_ind, phenotypes_cluster1] == "all"){
      pVals <- apply(norm_counts_subset,
                     1, function(mat){
                       x = mat[!grepl(phenotypes_cluster2, colnames(norm_counts_subset))]
                       y = mat[grepl(phenotypes_cluster2, colnames(norm_counts_subset))]
                       w = ks.test(x, y, alternative = "less")
                       c(w$p.value, w$statistic, median(x) - median(y))} )
    } else if(combinations[combinations_ind, phenotypes_cluster2] == "all"){
      pVals <- apply(norm_counts_subset,
                     1, function(mat){
                       x = mat[grepl(phenotypes_cluster1, colnames(norm_counts_subset))]
                       y = mat[!grepl(phenotypes_cluster1, colnames(norm_counts_subset))]
                       w = ks.test(x, y, alternative = "less")
                       c(w$p.value, w$statistic, median(x) - median(y))} )
    } else {
      pVals <- apply(norm_counts_subset,
                     1, function(mat){
                       x = mat[grepl(phenotypes_cluster1, colnames(norm_counts_subset))]
                       y = mat[grepl(phenotypes_cluster2, colnames(norm_counts_subset))]
                       w = ks.test(x, y, alternative = "less")
                       c(w$p.value, w$statistic, median(x) - median(y))} )
    }
  }
  data.table(phenotypes = combinations[combinations_ind, phenotypes],
             genes = colnames(pVals),
             pVals = pVals[1,],
             statistic = pVals[2,],
             diff_median = pVals[3,])
}

##' @rdname measureTrPhe
##' @name clusterCOMBS
##' @description \code{clusterCOMBS}: produces a data.table specifying all combinations of clusters for differential expression analysis.
##' @param clusters character vector, clusters (cell types) present in the data
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
