##' Summarise available perturbation data annotations using base::table
##' @rdname perturbTable
##' @name perturbTable
##' @author Vitalii Kleshchevnikov
##' @description \code{perturbTable} Cross Tabulates Connectivity map perturbation data annotations on columns of interest and subsets of perturbations
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param formula specify which 2 columns to call \code{\link[base]{table}} on: \code{pert_type ~ pert_time} or \code{type ~ pert_type} to tabulate single column
##' @param pert_ids data.table (optional) a way to select a subset of perturbation data annotations. data.table should be the output of \code{\link{geneName2PerturbAnno}} or \code{\link{geneID2PerturbAnno}} (should contain distil_id column, a unique perturbation ID)
##' @param pdata (optional) data.table, the output of \code{\link{openpData}}. By default \code{perturbTable} reads Cmap file (path provided in CMap_dirs), for multiple test it can be useful to supply data.table.
##' @return data.table containing the perturbation details from the Connectivity map project
##' @import data.table
##' @export perturbTable
##' @seealso \code{\link{openCellInfo}}, \code{\link{loadCMap}}, \code{\link{geneName2PerturbAnno}}
perturbTable = function(CMap_dirs, formula = pert_type ~ pert_time, pert_ids = NULL, pdata = NULL) {

  if(is.null(pdata)) pdata = openpData(CMap_dirs)

  if(length(formula) == 2){
    index1 = formula[[2]]
    if(is.null(pert_ids)) {
      table(pdata[,.(eval(index1))],
            dnn = c(as.character(index1)))
    } else {
      table(pdata[distil_id %in% pert_ids$distil_id,.(eval(index1))],
            dnn = c(as.character(index1)))
    }
  } else {
    index1 = formula[[2]]
    index2 = formula[[3]]
    if(is.null(pert_ids)) {
      table(pdata[,.(eval(index1), eval(index2))],
            dnn = c(as.character(index1), as.character(index2)))
    } else {
      table(pdata[distil_id %in% pert_ids$distil_id,.(eval(index1), eval(index2))],
            dnn = c(as.character(index1), as.character(index2)))
    }
  }
}
