##' Find regulators of gene sets using Connectivity Map perturbation data
##' @rdname regFromCMAP
##' @name regFromCMAP
##' @author Vitalii Kleshchevnikov
##' @description \code{regFromCMAP} identifies regulators of gene sets from Connectivity Map perturbation data using one-tailed Wilcox, KS or other (not yet implemented) tests
##' @param cmap object of class 'GCT' [package "cmapR"] with 7 slots containing z-score matrix, perturbation and feature details of the Connectivity map project. Returned by \code{\link{readCMAP}} or \code{\link{readCMAPsubset}}
##' @param gene_sets data.table containing gene set id and which genes are assigned to them
##' @param gene_set_id_col name of the column in \code{gene_sets} that stores gene set id
##' @param gene_id_col name of the column in \code{gene_sets} that stores gene id
##' @param method method for detecting differentially expressed genes. Currently only Wilcox and KS tests are implemented \link[stats]{wilcox.test}, \link[stats]{ks.test}.
##' @param GSEA_weighting Weighting of ranking/correlations, see \link[gsEasy]{gset}, https://cran.r-project.org/web/packages/gsEasy/vignettes/gsEasy-guide.html or Subramanian et. al 2005.
##' @param cutoff FDR-corrected p-value cutoff
##' @param pval_corr_method multiple hypothesis p-value correction method. Details: \link[stats]{p.adjust} - method.
##' @param renormalise renormalise z-scores for a Connectivity map subset being analysed. Z-score = (X-\link[stats]{median}(X)) / (\link[stats]{mad}(X)*1.4826)
##' @param n_cores number of cores to be used in parallel processing (over combinations of clusters). More details: \link[parallel]{parLapply}, \link[parallel]{makeCluster}, \link[parallel]{detectCores}
##' @return data.table containing gene set id, which genes may regulate this set, test statistic and difference in medians between each gene set and all other genes
##' @import data.table
##' @import parallel
##' @export regFromCMAP
##' @examples
##' library(regNETcmap)
##' library(R.utils)
##' # read all knockdown perturbations (3.4GB)
##' CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
##' # first, let's find cell lines and measurement times
##' CMAPsub = readCMAPsubset(is_touchstone = T, pert_types = c("trt_sh.cgs"), pert_times = c("all"), cell_ids = c("all"), CMap_files)
##' pdata = as.data.table(CMAPsub@cdesc)
##' pdata$cell_id = as.character(pdata$cell_id)
##' pdata$pert_time = as.character(pdata$pert_time)
##' pdata$pert_type = as.character(pdata$pert_type)
##' perturbTable(formula = cell_id ~ pert_time, pdata = pdata)
##' # Pick one cell line and measurement time
##' CMAPsub = readCMAPsubset(is_touchstone = T, pert_types = c("trt_sh.cgs"), pert_times = c("96"), cell_ids = c("A375"), CMap_files)
##' # Read gene sets
##' trPhe_pair_map_file = "../regulatory_networks_by_cmap/results/phenotypes_genes_pvals_pairwise_mapped"
##' trPhe_pair_map_file.zip = paste0(trPhe_pair_map_file, ".gz")
##' gunzip(trPhe_pair_map_file.zip, overwrite = T, remove = F)
##' pVals_human = fread(trPhe_pair_map_file, sep = "\t", stringsAsFactors = F)
##' unlink(trPhe_pair_map_file)
##' # Select genes based on corrected p-value and difference in medians
##' pVals_human = pVals_human[pVals_corr < 0.05 & diff_median > 0]
##' pVals_reg = regFromCMAP(cmap = CMAPsub,
##'     gene_sets = pVals_human,
##'     gene_set_id_col = "phenotypes", gene_id_col = "entrezgene",
##'     method = "ks", cutoff = 1, pval_corr_method = "fdr",
##'     n_cores = detectCores() - 1)
##' qplot(x = pVals_reg$diff_median, y = -log10(pVals_reg$pVals), geom = "bin2d", bins = 150) + theme_light()
regFromCMAP = function(cmap, gene_sets, gene_set_id_col = "phenotypes", gene_id_col = "entrezgene", method =  c("ks", "wilcox", "GSEA")[1], GSEA_weighting = 1, cutoff = 1, pval_corr_method = "fdr", renormalise = F, n_cores = detectCores() - 1){
  gene_sets = unique(gene_sets[, c(gene_set_id_col, gene_id_col), with = F])
  setnames(gene_sets, c(gene_set_id_col, gene_id_col), c("gene_set_id", "gene_id"))
  setorder(gene_sets, gene_set_id)
  # remove not mapped genes
  gene_sets = gene_sets[gene_id != "" | !is.na(gene_id) | gene_id != "NA"]
  # extract data matrix
  cmap_mat = cmap@mat
  # renormalise:
  if(renormalise) cmap_mat = t(apply(cmap_mat, 1, function(X) (X-median(X)) / (mad(X)*1.4826)))
  # extract column description
  cmap_cdesc = as.data.table(cmap@cdesc)

  # determine if shRNA or overexpression are being tested, modify negative vs positive
  if(length(unique(cmap_cdesc$pert_type)) != 1) stop("regFromCMAP: only single perturbation type (pert_type) is allowed")
  if(unique(cmap_cdesc$pert_type)[1] == "trt_oe") {
    alternative_neg = "less"
    alternative_pos = "greater"
  } else { #unique(cmap_cdesc$pert_type) == "trt_sh.cgs" or anything else
    alternative_neg = "greater"
    alternative_pos = "less"
  }
  # determine if ks or wilcox test is being used, modify negative vs positive
  if(method == "wilcox"){
    alternative_neg_temp = alternative_pos
    alternative_pos_temp = alternative_neg
    alternative_neg = alternative_neg_temp
    alternative_pos = alternative_pos_temp
  }

  # set up parallel processing
  # create cluster
  cl <- makeCluster(n_cores)
  # get library support needed to run the code
  clusterEvalQ(cl, {library(regNETcmap)})
  # put objects in place that might be needed for the code
  clusterExport(cl, c("cmap_mat", "gene_sets", "method", "alternative_neg", "alternative_pos", "GSEA_weighting"), envir=environment())

  pVals_neg = parLapply(cl, X = unique(gene_sets$gene_set_id),
                        fun = function(name) {
                          regFromCMAPSingle(cmap_mat = cmap_mat,
                                            gene_sets = gene_sets,
                                            gene_set_name = name,
                                            method = method,
                                            GSEA_weighting = GSEA_weighting,
                                            alternative = alternative_neg)
                        })
  pVals_pos = parLapply(cl, X = unique(gene_sets$gene_set_id),
                        fun = function(name) {
                          regFromCMAPSingle(cmap_mat = cmap_mat,
                                            gene_sets = gene_sets,
                                            gene_set_name = name,
                                            method = method,
                                            GSEA_weighting = GSEA_weighting,
                                            alternative = alternative_pos)
                        })
  # stop cluster
  stopCluster(cl)

  # combine results
  pVals_neg = Reduce(rbind, pVals_neg)
  pVals_neg$sign = "negative"
  pVals_pos = Reduce(rbind, pVals_pos)
  pVals_pos$sign = "positive"
  pVals = rbind(pVals_neg, pVals_pos)
  # map signatures to genes
  cmap_cdesc_temp = cmap_cdesc[,.(sig_id, pert_iname, pert_type, cell_id, pert_time, pert_idose)]
  pVals = merge(x = pVals, y = cmap_cdesc_temp,
                by = "sig_id",
                all.x = T, all.y = F, allow.cartesian = TRUE)
  # add total number of genes in a set
  genes_in_set = gene_sets[, .(genes_in_set = uniqueN(gene_id)), by = "gene_set_id"]
  pVals = merge(x = pVals, y = genes_in_set,
                by.x = "gene_set_id", by.y = "gene_set_id",
                all.x = T, all.y = F, allow.cartesian = TRUE)
  # multiple testing correction
  pVals[, pVals_corr := p.adjust(pVals, method = pval_corr_method)]
  pVals[, passed := pVals_corr < cutoff]
  pVals
}

##' @rdname regFromCMAP
##' @name regFromCMAPSingle
##' @description \code{regFromCMAPSingle}: identifies regulators of a single gene set from Connectivity Map perturbation data using one-tailed Wilcox, KS or other (not yet implemented) tests
##' @param cmap_mat Connectivity map data matrix (\code{cmap@mat})
##' @param gene_set_name name of one gene set in gene_sets
##' @param alternative alternative hypothesis for statistical tests, \link[stats]{wilcox.test}, \link[stats]{ks.test}
##' @return \code{regFromCMAPSingle}: data.table containing gene set id (single), which genes may regulate this set, test statistic and difference in medians between each gene set and all other genes
##' @import data.table
##' @export regFromCMAPSingle
regFromCMAPSingle = function(cmap_mat, gene_sets,
                             gene_set_name,
                             method = c("ks","wilcox", "GSEA")[1],
                             GSEA_weighting = 1,
                             alternative = c("less", "greater")[1]){

  genes = gene_sets[gene_set_id == gene_set_name, gene_id]
  genes_ind = rownames(cmap_mat) %in% genes

  if(method == "wilcox"){
    pVals <- apply(cmap_mat,
                   2, function(mat){
                     x = mat[!genes_ind]
                     y = mat[genes_ind]
                     w = wilcox.test(x, y, alternative = alternative)
                     c(w$p.value, w$statistic, median(y) - median(x))} )
  }
  if(method == "ks"){
    pVals <- apply(cmap_mat,
                   2, function(mat){
                     x = mat[!genes_ind]
                     y = mat[genes_ind]
                     w = ks.test(x, y, alternative = alternative)
                     c(w$p.value, w$statistic, median(y) - median(x))} )
  }
  if(method == "GSEA"){
    if(alternative == "less"){
      x = proc.time()
      pVals <- apply(cmap_mat,
                     2, function(mat){
                       p.val = gsEasy::gset(S = which(genes_ind), N = NULL,
                                            r = -mat, p = GSEA_weighting, min_its = 1000,
                                            max_its = 1e+05,
                                            significance_threshold = 1, log_dismiss = -10,
                                            raw_score = FALSE)
                       score = gsEasy::gset(S = which(genes_ind), N = NULL,
                                            r = -mat, p = GSEA_weighting, min_its = 1000,
                                            max_its = 1e+05,
                                            significance_threshold = 1, log_dismiss = -10,
                                            raw_score = TRUE)
                       c(p.val, score, median(mat[genes_ind]) - median(mat[!genes_ind]))} )
      proc.time() - x
    } else if (alternative == "greater"){
      pVals <- apply(cmap_mat,
                     2, function(mat){
                       p.val = gsEasy::gset(S = which(genes_ind), N = NULL,
                                            r = mat, p = GSEA_weighting, min_its = 1000,
                                            max_its = 1e+05,
                                            significance_threshold = 1, log_dismiss = -10,
                                            raw_score = FALSE)
                       score = gsEasy::gset(S = which(genes_ind), N = NULL,
                                            r = mat, p = GSEA_weighting, min_its = 1000,
                                            max_its = 1e+05,
                                            significance_threshold = 1, log_dismiss = -10,
                                            raw_score = TRUE)
                       c(p.val, score, median(mat[genes_ind]) - median(mat[!genes_ind]))} )
    }

  }

  data.table(gene_set_id = gene_set_name,
             sig_id = colnames(pVals),
             pVals = pVals[1,],
             statistic = pVals[2,],
             diff_median = pVals[3,],
             genes_in_cmap = sum(genes_ind))
}

##' @rdname regFromCMAP
##' @name regFromCMAPmcell
##' @description \code{regFromCMAPmcell}: identifies regulators of gene sets from multiple cell lines Connectivity Map perturbation data using one-tailed Wilcox, KS or other (not yet implemented) tests
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param cell_ids a character vector of cell line names, query for available cell line names using \code{perturbTable(CMap_files, ~ cell_id)}
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}
##' @param pert_times a character vector of perturbation times, query for available perturbation times using \code{perturbTable(CMap_files, ~ pert_time)}
##' @param is_touchstone logical, select only the perturbation data that is a part of the touchstone dataset (genes profiled across all 9 core cell lines)? Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/
##' @param min_cell_lines keep a gene=>gene-set interaction if at least \code{min_cell_lines} support it (inclusive)
##' @param max_cell_lines keep a gene=>gene-set interaction if at most \code{max_cell_lines} support it (inclusive)
##' @param gene_names a character vector of HUGO gene names. Use \code{"all"} to select all perturbations.
##' @param keep_one_oe keep only one overexpression experiment per gene and condition. Applicable only when \code{pert_types = "trt_oe"}. Perturbations where sig_id matches pert_id are retained ("one"). To invert the selection use "other". To select all use "all"
##' @return \code{regFromCMAPmcell}: data.table containing gene set id, which genes may regulate these sets and the cell line of origin, test statistic and difference in medians between each gene set and all other genes
##' @import data.table
##' @export regFromCMAPmcell
regFromCMAPmcell = function(CMap_files, cell_ids = c("A375", "A549", "HA1E",
                                                     "HCC515", "HEPG2", "HT29",
                                                     "MCF7", "PC3"),
                            pert_types = c("trt_sh.cgs", "trt_oe", "trt_xpr"),#"ctl_untrt.cns", "ctl_vector.cns", "ctl_vehicle.cns"),
                            pert_times = c("96"),
                            is_touchstone = T, # F (ctl is not touchstone)
                            min_cell_lines = 1, max_cell_lines = 9,
                            gene_sets,
                            gene_set_id_col = "phenotypes",
                            gene_id_col = "entrezgene",
                            method =  c("ks", "wilcox", "GSEA")[1],
                            GSEA_weighting = 1,
                            keep_one_oe = c("one", "other", "all")[1],
                            cutoff = 0.05, pval_corr_method = "fdr",
                            renormalise = F,
                            n_cores = detectCores() - 1,
                            gene_names = c("all", unique(gene_sets$hgnc_symbol)[3:22])[1]) {
  if(!min_cell_lines <= max_cell_lines) stop("regFromCMAPmcell: min_cell_lines should be smaller or equal to max_cell_lines")
  if(length(pert_times) != 1) stop("regFromCMAPmcell: only single perturbation time (pert_times) is allowed")
  types = perturbTable(CMap_files, ~ pert_type)
  pert_types = pert_types[pert_types %in% names(types)]
  lines = perturbTable(CMap_files, ~ cell_id)
  cell_ids = cell_ids[cell_ids %in% names(lines)]

  pVals_cell = lapply(cell_ids, function(cell_line, pert_types){
    pVals_type = lapply(pert_types, function(pert_type, cell_line){
      CMAPsub = readCMAPsubset(is_touchstone = is_touchstone, pert_types = pert_type,
                               pert_times = pert_times, cell_ids = cell_line,
                               CMap_files, gene_names = gene_names, keep_one_oe = keep_one_oe)
      regFromCMAP(cmap = CMAPsub,
                  gene_sets = gene_sets,
                  gene_set_id_col = gene_set_id_col, gene_id_col = gene_id_col,
                  method = method, GSEA_weighting = GSEA_weighting,
                  cutoff = 1, pval_corr_method = pval_corr_method,
                  renormalise = renormalise,
                  n_cores = n_cores)
    }, cell_line)
    # combine results
    Reduce(rbind, pVals_type)
  }, pert_types)
  # combine results
  pVals_cell = Reduce(rbind, pVals_cell)
  # multiple testing correction
  pVals_cell[, pVals_corr := p.adjust(pVals, method = pval_corr_method)]
  pVals_cell[, passed := pVals_corr < cutoff]
  pVals_cell_thresh = copy(pVals_cell[passed == T])

  # support from multiple cell lines
  pVals_cell_thresh[, n_cell_id := uniqueN(cell_id), by = .(gene_set_id, pert_iname, pert_type)]
  pVals_cell_thresh = pVals_cell_thresh[n_cell_id >= min_cell_lines & n_cell_id <= max_cell_lines]

  list(pVals_cell = pVals_cell,
       pVals_cell_thresh = pVals_cell_thresh)
}
