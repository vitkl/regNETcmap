##' Read a subsets of perturbations (Connectivity Map data, all pertubages, different conditions)
##' @rdname readCMAPsubset
##' @name readCMAPsubset
##' @author Vitalii Kleshchevnikov
##' @description \code{readCMAPsubset} uses \code{\link{geneName2PerturbAnno}} to select all perturbations in a set of cell lines, perturbation types (compound, shRNA, overexpression, e.g.) and time and then uses \code{readCMAP} to read Connectivity Map data given these conditions.
##' @param is_touchstone logical, select only the perturbation data that is a part of the touchstone dataset (genes profiled across all 9 core cell lines)? Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}
##' @param pert_times a character vector of perturbation times, query for available perturbation times using \code{perturbTable(CMap_files, ~ pert_time)}
##' @param cell_ids a character vector of cell line names, query for available cell line names using \code{perturbTable(CMap_files, ~ cell_id)}
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return object of class 'GCT' [package "cmapR"] with 7 slots containing z-score matrix, perturbation and feature details of the Connectivity map project
##' @export readCMAPsubset
##' @seealso \code{\link{openCellInfo}}, \code{\link{loadCMap}}, \code{\link{perturbTable}}
##' @examples
##' # read all knockdown perturbations (3.4GB)
##' CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
##' allIDs = geneName2PerturbAnno(gene_names = "all", CMap_files = CMap_files, is_touchstone = T, pert_types = "trt_sh.cgs")
##' perturbTable(pdata = allIDs)
##' CMAP = readCMAP(PerturbAnno = allIDs, CMap_files = CMap_files)
##' # fast way to read the same data
##' CMAPsub = readCMAPsubset(is_touchstone = T, pert_types = c("trt_sh.cgs"), pert_times = c("all"), cell_ids = c("all"), CMap_files)
##' perturbTable(pdata = as.data.table(CMAPsub@cdesc))
readCMAPsubset = function(is_touchstone = c("all", T, F)[2],
                          pert_types = c("trt_sh.cgs", "trt_oe"),
                          pert_times = c("all"),
                          cell_ids = c("all"), CMap_files) {
  allIDs = geneName2PerturbAnno(gene_names = "all", CMap_files,
                                is_touchstone = is_touchstone, pert_types = pert_types,
                                pert_times = pert_times, cell_ids = cell_ids)
  readCMAP(PerturbAnno = allIDs, CMap_files = CMap_files)
}