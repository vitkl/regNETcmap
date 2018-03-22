##' Summarise TF targets enrichment when TF perturbed in Connectivity Map
##' @rdname TFtargetsINcmap
##' @name TFtargetsINcmap
##' @author Vitalii Kleshchevnikov
##' @description \code{TFtargetsINcmap}
##' @param regulons data.table containing TF regulons (TF and target columns)
##' @param alternative indicates alternative hypothesis for \link[stats]{ks.test}. "trt_oe": negative_regulator = "less", positive_regulator = "greater"; "trt_sh.cgs" or anything else: negative_regulator = "greater"; positive_regulator = "less"
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}
##' @param cell_ids a character vector of cell line names, query for available cell line names using \code{perturbTable(CMap_files, ~ cell_id)}
##' @param gene_names a character vector of HUGO gene names. Use \code{"all"} to select all perturbations.
##' @param pert_times a character vector of perturbation times, query for available perturbation times using \code{perturbTable(CMap_files, ~ pert_time)}
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param keep_one_oe keep only one overexpression experiment per gene and condition. Applicable only when \code{pert_types = "trt_oe"}. Perturbations where sig_id matches pert_id are retained ("one"). To invert the selection use "other". To select all use "all"
##' @param landmark_only look only at landmark genes. Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit
##' @param clustermq_memory When using clustermq: memory requested for each job
##' @param clustermq_job_size When using clustermq: The number of function calls per job
##' @param min_targets_in_cmap minimal number of TF targets in cmap
##' @param force_gene_names look at all \code{gene_names} but less cell lines
##' @return list containing data.table summary of TF perturbation effects on TF targets, object of class 'GCT' (CMap data), summary of genes vs cell lines available in CMap and returned (gene_cell_counts)
##' @export TFtargetsINcmap
##' @import data.table
##' @import clustermq
##' @import gsEasy
##' @seealso \code{\link{regFromCMAP}}, \code{\link{completeCellGeneCombs}}, \code{\link{showLargestRegulons}}
TFtargetsINcmap = function(regulons, alternative = "less",
                           pert_types = c("trt_oe"),
                           cell_ids = c("A375", "A549", "HA1E","HCC515",
                                        "HEPG2", "HT29","MCF7", "PC3"),
                           gene_names = c("NFKB1", "E2F1", "MYC", "TP53",
                                          "YY1", "JUN", "STAT1", "STAT3"),
                           pert_times = c("96"),
                           CMap_files, keep_one_oe = c("one", "other", "all")[1],
                           landmark_only = F,
                           clustermq_memory = 2000, clustermq_job_size = 1,
                           min_targets_in_cmap = 3, force_gene_names = F
){


  cmap = readCMAPsubset(is_touchstone = "all",
                        pert_types = pert_types, #trt_sh.cgs trt_oe
                        pert_times = pert_times, cell_ids = cell_ids,
                        gene_names = gene_names, CMap_files = CMap_files,
                        keep_one_oe = keep_one_oe,
                        landmark_only = landmark_only)

  gene_cell_counts = completeCellGeneCombs(cmap, force_gene_names)

  cell_ids = cell_ids[cell_ids %in% colnames(gene_cell_counts$gene_cell_counts)]
  gene_names = gene_names[gene_names %in% rownames(gene_cell_counts$gene_cell_counts)]

  # look only at regulons where target genes were measured in cmap
  regulons = regulons[target_entrezgene %in% rownames(cmap@mat)]
  regulons[, targets_in_cmap := uniqueN(target_entrezgene), by = TF]
  regulons = regulons[targets_in_cmap >= min_targets_in_cmap]
  gene_names = gene_names[gene_names %in% regulons$TF]

  # produce data for plotting
  res = Q(function(cell_line){
    library(regNETcmap)
    library(gsEasy)
    res = lapply(gene_names, function(TF_sh){
      res = lapply(gene_names, function(TF_measured){
        mat = cmap@mat
        target_ind = rownames(mat) %in% regulons[TF %in% TF_measured, target_entrezgene]
        cell_line_ind = cmap@cdesc$cell_id %in% cell_line
        TF_sh_ind = cmap@cdesc$pert_iname %in% TF_sh
        x = mat[!target_ind, cell_line_ind & TF_sh_ind]
        y = mat[target_ind, cell_line_ind & TF_sh_ind]
        w = ks.test(x, y, alternative = alternative)
        if(alternative == "less") GSEA_mat = -mat else GSEA_mat = mat
        GSEA_pval1 = gsEasy::gset(S = which(target_ind), N = NULL,
                                  r = GSEA_mat[, cell_line_ind & TF_sh_ind], p = 1, min_its = 1000,
                                  max_its = 1e+05,
                                  significance_threshold = 1, log_dismiss = -10,
                                  raw_score = FALSE)
        #GSEA_score1 = gsEasy::gset(S = which(target_ind), N = NULL,
        #                           r = GSEA_mat, p = 1, min_its = 1000,
        #                           max_its = 1e+05,
        #                           significance_threshold = 1, log_dismiss = -10,
        #                           raw_score = TRUE)
        GSEA_pval10 = gsEasy::gset(S = which(target_ind), N = NULL,
                                   r = GSEA_mat[, cell_line_ind & TF_sh_ind], p = 10, min_its = 1000,
                                   max_its = 1e+05,
                                   significance_threshold = 1, log_dismiss = -10,
                                   raw_score = FALSE)
        #GSEA_score10 = gsEasy::gset(S = which(target_ind), N = NULL,
        #                            r = GSEA_mat, p = 10, min_its = 1000,
        #                            max_its = 1e+05,
        #                            significance_threshold = 1, log_dismiss = -10,
        #                            raw_score = TRUE)
        size = sum(target_ind)
        data.table(target = ifelse(target_ind,"TF_targets", "other_genes"), TF_sh_lab = paste0(TF_sh," shRNA"),
                   cell_ids = cell_line, TF_measured_lab = paste0(TF_measured, " regulon\ntargets:", size),
                   TF_measured = TF_measured, TF_sh = TF_sh,
                   gene_exp = mat[, cell_line_ind & TF_sh_ind],
                   pval = w$p.value, statistic = w$statistic,
                   GSEA_pval1 = GSEA_pval1, GSEA_pval10 = GSEA_pval10,
                   #GSEA_score1 = GSEA_score1, GSEA_score10 = GSEA_score10,
                   size = regulons[TF %in% TF_measured, unique(size)]
        )
      })
      Reduce(rbind, res)
    })
    Reduce(rbind, res)
  }, cell_ids,
  export = list(regulons = regulons, gene_names = gene_names, cmap = cmap, alternative = alternative),
  job_size = clustermq_job_size, memory = clustermq_memory)
  res = Reduce(rbind, res)
  res[, pvals := paste0("ks: ", signif(pval, 3), "\n",
                        signif(GSEA_pval1, 3), "\n",
                        signif(GSEA_pval10, 3)
                        )]
  list(res = res, cmap = cmap, gene_cell_counts = gene_cell_counts, regulons = regulons)
}

##' Plot TF perturbation effects (CMap) on TF targets
##' @rdname plotTFtargetsINcmap
##' @name plotTFtargetsINcmap
##' @param res data.table containing summary of TF perturbation effects on TF targets returned by \code{\link{TFtargetsINcmap}}
##' @param TFsels which TF to plot, defaults to all in \code{res}
##' @import data.table
##' @export plotTFtargetsINcmap
plotTFtargetsINcmap = function(res, TFsels = NULL,
                               title = "TF that positively regulate their targets \n(shRNA induces shift to the left)",
                               lab_text_size = 3, lab_text_y = 1.5, lab_text_x = 5, strip.text_size = 8) {
  if(is.null(TFsels)) TFsels = unique(c(res$TF_sh, res$TF_measured))
  plots = list()
  for (cellline in unique(res$cell_ids)) {
    p = ggplot(res[cell_ids == cellline & TF_sh %in% TFsels & TF_measured %in% TFsels],
               aes(x = gene_exp, color = target)) +
      geom_density() + facet_grid(TF_measured_lab~TF_sh_lab) +
      geom_text(y = lab_text_y, x = lab_text_x, aes(label = pvals), size = lab_text_size) +
      theme(strip.text.y = element_text(angle = 0, size = strip.text_size),
            strip.text.x = element_text(size = strip.text_size)) +
      ggtitle(title, subtitle = cellline) + xlab("z-score")
    plots[[which(unique(res$cell_ids) %in% cellline)]] = p
  }
  plots
}

##' Show largest TF regulons
##' @rdname showLargestRegulons
##' @name showLargestRegulons
##' @param regulons data.table containing TF regulons (TF and target columns)
##' @param N number of largest regulons to display
##' @import data.table
##' @export showLargestRegulons
showLargestRegulons = function(regulons, N = 40) {
  regulons = copy(regulons)
  regulons[, size := uniqueN(target), by = TF]
  regulons = regulons[order(size, decreasing = T)]
  list(regulons = regulons, summary = unique(regulons[order(size, decreasing = T),.(TF, size)])[1:N])
}

##' Find cell lines and perturbations (Connectivity map) where all pert_iname are present in all cell_id
##' @rdname completeCellGeneCombs
##' @name completeCellGeneCombs
##' @param cmap object of class 'GCT' [package "cmapR"] with 7 slots containing z-score matrix, perturbation and feature details of the Connectivity map project. Returned by \code{\link{readCMAP}} or \code{\link{readCMAPsubset}}
##' @param force_gene_names look at all \code{gene_names} but less cell lines
##' @import data.table
##' @export completeCellGeneCombs
completeCellGeneCombs = function(cmap, force_gene_names = F) {
  gene_cell = as.data.table(cmap@cdesc)[, .(pert_iname, cell_id)]
  gene_cell$pert_iname = as.character(gene_cell$pert_iname)
  gene_cell$cell_id = as.character(gene_cell$cell_id)
  gene_cell_counts = table(gene_cell)
  if(ncol(gene_cell_counts) > 1 & nrow(gene_cell_counts) > 1){
    gene_zero_counts = rowMeans(gene_cell_counts == 0)
    cell_zero_counts = colMeans(gene_cell_counts == 0)
    if(max(gene_zero_counts) < max(cell_zero_counts) | force_gene_names){
      gene_cell_counts = gene_cell_counts[,c(!as.logical(cell_zero_counts))]
      gene_zero_counts = rowMeans(gene_cell_counts == 0)
      cell_zero_counts = colMeans(gene_cell_counts == 0)
      gene_cell_counts = gene_cell_counts[!as.logical(gene_zero_counts),]
      gene_cell_counts = gene_cell_counts[,!as.logical(cell_zero_counts)]
    } else {
      gene_cell_counts = gene_cell_counts[!as.logical(gene_zero_counts),]
      gene_zero_counts = rowMeans(gene_cell_counts == 0)
      cell_zero_counts = colMeans(gene_cell_counts == 0)
      gene_cell_counts = gene_cell_counts[!as.logical(gene_zero_counts),]
      gene_cell_counts = gene_cell_counts[,!as.logical(cell_zero_counts)]
    }
  }
  list(gene_cell_counts = gene_cell_counts, gene_cell_counts_requested = table(gene_cell))
}
