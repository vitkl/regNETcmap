% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regFromCMAP.R
\name{regFromCMAP}
\alias{regFromCMAP}
\alias{regFromCMAPSingle}
\alias{regFromCMAPmcell}
\title{Find regulators of gene sets using Connectivity Map perturbation data}
\usage{
regFromCMAP(cmap, gene_sets, gene_set_id_col = "phenotypes",
  gene_id_col = "entrezgene", method = c("ks", "wilcox", "GSEA")[1],
  GSEA_weighting = 1, cutoff = 1, pval_corr_method = "fdr",
  renormalise = F, n_cores = detectCores() - 1, clustermq = F,
  clustermq_seed = 128965, clustermq_memory = 2000,
  clustermq_job_size = 10)

regFromCMAPSingle(cmap_mat, gene_sets, gene_set_name, method = c("ks",
  "wilcox", "GSEA")[1], GSEA_weighting = 1, alternative = c("less",
  "greater")[1])

regFromCMAPmcell(CMap_files, cell_ids = c("A375", "A549", "HA1E",
  "HCC515", "HEPG2", "HT29", "MCF7", "PC3"), pert_types = c("trt_sh.cgs",
  "trt_oe", "trt_xpr"), pert_itimes = c("96 h"), is_touchstone = T,
  min_cell_lines = 1, max_cell_lines = 9, gene_sets,
  gene_set_id_col = "phenotypes", gene_id_col = "entrezgene",
  method = c("ks", "wilcox", "GSEA")[1], GSEA_weighting = 1,
  keep_one_oe = c("one", "other", "all")[1], cutoff = 0.05,
  pval_corr_method = "fdr", renormalise = F, n_cores = detectCores()
  - 1, gene_names = c("all", unique(gene_sets$hgnc_symbol)[3:22])[1],
  landmark_only = F)
}
\arguments{
\item{cmap}{object of class 'GCT' [package "cmapR"] with 7 slots containing z-score matrix, perturbation and feature details of the Connectivity map project. Returned by \code{\link{readCMAP}} or \code{\link{readCMAPsubset}}}

\item{gene_sets}{data.table containing gene set id and which genes are assigned to them}

\item{gene_set_id_col}{name of the column in \code{gene_sets} that stores gene set id}

\item{gene_id_col}{name of the column in \code{gene_sets} that stores gene id}

\item{method}{method for detecting differentially expressed genes. Currently only Wilcox and KS tests are implemented \link[stats]{wilcox.test}, \link[stats]{ks.test}.}

\item{GSEA_weighting}{Weighting of ranking/correlations, see \link[gsEasy]{gset}, https://cran.r-project.org/web/packages/gsEasy/vignettes/gsEasy-guide.html or Subramanian et. al 2005.}

\item{cutoff}{FDR-corrected p-value cutoff}

\item{pval_corr_method}{multiple hypothesis p-value correction method. Details: \link[stats]{p.adjust} - method.}

\item{renormalise}{renormalise z-scores for a Connectivity map subset being analysed. Z-score = (X-\link[stats]{median}(X)) / (\link[stats]{mad}(X)*1.4826)}

\item{n_cores}{number of cores to be used in parallel processing (over combinations of clusters). More details: \link[parallel]{parLapply}, \link[parallel]{makeCluster}, \link[parallel]{detectCores}}

\item{clustermq}{Use clustermq LSF job scheduler (TRUE) as an alternative to parLapply (FALSE). Details: \link[clustermq]{Q}}

\item{clustermq_seed}{When using clustermq: Seed for random number generation.}

\item{clustermq_memory}{When using clustermq: memory requested for each job}

\item{clustermq_job_size}{When using clustermq: The number of function calls per job}

\item{cmap_mat}{Connectivity map data matrix (\code{cmap@mat})}

\item{gene_set_name}{name of one gene set in gene_sets}

\item{alternative}{alternative hypothesis for statistical tests, \link[stats]{wilcox.test}, \link[stats]{ks.test}}

\item{CMap_files}{a list of directories and urls produced by \code{\link{loadCMap}}}

\item{cell_ids}{a character vector of cell line names, query for available cell line names using \code{perturbTable(CMap_files, ~ cell_id)}}

\item{pert_types}{a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}}

\item{pert_itimes}{a character vector of perturbation times, query for available perturbation times using \code{perturbTable(CMap_files, ~ pert_itime)}}

\item{is_touchstone}{logical, select only the perturbation data that is a part of the touchstone dataset (genes profiled across all 9 core cell lines)? Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/}

\item{min_cell_lines}{keep a gene=>gene-set interaction if at least \code{min_cell_lines} support it (inclusive)}

\item{max_cell_lines}{keep a gene=>gene-set interaction if at most \code{max_cell_lines} support it (inclusive)}

\item{keep_one_oe}{keep only one overexpression experiment per gene and condition. Applicable only when \code{pert_types = "trt_oe"}. Perturbations where sig_id matches pert_id are retained ("one"). To invert the selection use "other". To select all use "all"}

\item{gene_names}{a character vector of HUGO gene names. Use \code{"all"} to select all perturbations.}

\item{landmark_only}{look only at landmark genes. Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit}
}
\value{
data.table containing gene set id, which genes may regulate this set, test statistic and difference in medians between each gene set and all other genes

\code{regFromCMAPSingle}: data.table containing gene set id (single), which genes may regulate this set, test statistic and difference in medians between each gene set and all other genes

\code{regFromCMAPmcell}: data.table containing gene set id, which genes may regulate these sets and the cell line of origin, test statistic and difference in medians between each gene set and all other genes
}
\description{
\code{regFromCMAP} identifies regulators of gene sets from Connectivity Map perturbation data using one-tailed Wilcox, KS or other (not yet implemented) tests

\code{regFromCMAPSingle}: identifies regulators of a single gene set from Connectivity Map perturbation data using one-tailed Wilcox, KS or other (not yet implemented) tests

\code{regFromCMAPmcell}: identifies regulators of gene sets from multiple cell lines Connectivity Map perturbation data using one-tailed Wilcox, KS or other (not yet implemented) tests
}
\examples{
library(regNETcmap)
library(R.utils)
# read all knockdown perturbations (3.4GB)
CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
# first, let's find cell lines and measurement times
CMAPsub = readCMAPsubset(is_touchstone = T, pert_types = c("trt_sh.cgs"), pert_itimes = c("all"), cell_ids = c("all"), CMap_files)
pdata = as.data.table(CMAPsub@cdesc)
pdata$cell_id = as.character(pdata$cell_id)
pdata$pert_itime = as.character(pdata$pert_itime)
pdata$pert_type = as.character(pdata$pert_type)
perturbTable(formula = cell_id ~ pert_itime, pdata = pdata)
# Pick one cell line and measurement time
CMAPsub = readCMAPsubset(is_touchstone = T, pert_types = c("trt_sh.cgs"), pert_itimes = c("96 h"), cell_ids = c("A375"), CMap_files)
# Read gene sets
trPhe_pair_map_file = "../regulatory_networks_by_cmap/results/phenotypes_genes_pvals_pairwise_mapped"
trPhe_pair_map_file.zip = paste0(trPhe_pair_map_file, ".gz")
gunzip(trPhe_pair_map_file.zip, overwrite = T, remove = F)
pVals_human = fread(trPhe_pair_map_file, sep = "\\t", stringsAsFactors = F)
unlink(trPhe_pair_map_file)
# Select genes based on corrected p-value and difference in medians
pVals_human = pVals_human[pVals_corr < 0.05 & diff_median > 0]
pVals_reg = regFromCMAP(cmap = CMAPsub,
    gene_sets = pVals_human,
    gene_set_id_col = "phenotypes", gene_id_col = "entrezgene",
    method = "ks", cutoff = 1, pval_corr_method = "fdr",
    n_cores = detectCores() - 1)
qplot(x = pVals_reg$diff_median, y = -log10(pVals_reg$pVals), geom = "bin2d", bins = 150) + theme_light()
}
\author{
Vitalii Kleshchevnikov
}
