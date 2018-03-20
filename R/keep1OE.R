##' Keep only 1 overexpression experiment per gene and condition
##' @rdname geneName2PerturbAnno
##' @name keep1OE
##' @author Vitalii Kleshchevnikov
##' @description \code{keep1OE} Keep only 1 overexpression experiment per gene and condition. Designed to be applicable to \code{pert_types = "trt_oe"}.
##' @param PerturbAnno data.table containing the perturbation details from the Connectivity map project produced by \code{\link{geneName2PerturbAnno}}
##' @param keep_one_oe keep only one overexpression experiment per gene and condition. Applicable only when \code{pert_types = "trt_oe"}. Perturbations where sig_id matches pert_id are retained ("one"). To invert the selection use "other". To select all use "all"
##' @param pert_types find duplicates in this perturbation type
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return data.table containing the perturbation details from the Connectivity map project. Only one overexpression experiment per gene and condition is retained.
##' @import data.table
##' @export keep1OE
keep1OE = function(PerturbAnno, keep_one_oe = c("one", "other", "all")[1], pert_types = "trt_oe", CMap_files) {
  duplicated = PerturbAnno[pert_type == pert_types[1]]
  duplicated[, n_per_iname := uniqueN(sig_id), by = .(pert_iname, pert_type, cell_id, pert_time)]
  duplicated = duplicated[n_per_iname > 1]
  if(nrow(duplicated) >= 1){
    duplicated_IDs = geneName2PerturbAnno(gene_names = unique(duplicated$pert_iname), CMap_files = CMap_files,
                                          is_touchstone = as.logical(unique(duplicated$is_touchstone)),
                                          pert_types = pert_types[1],
                                          pert_times = unique(as.character(duplicated$pert_time)),
                                          cell_ids = unique(as.character(duplicated$cell_id)))
    duplicated_IDs = duplicated_IDs[order(pert_iname)]
    duplicated_IDs[, pert_id_in_sig_id := grep(pert_id, sig_id, ignore.case = T, value = T), by = pert_id]
    if(keep_one_oe == "one"){
      duplicated_IDs = duplicated_IDs[sig_id == pert_id_in_sig_id]
    } else if(keep_one_oe == "other") {
      duplicated_IDs = duplicated_IDs[sig_id != pert_id_in_sig_id]
    } else message("keep_one_oe argument can be only \"one\", \"other\", \"all\". No perturbation filtered")
    duplicated_IDs$pert_id_in_sig_id = NULL
    # add duplicated to all
    PerturbAnno = PerturbAnno[!(pert_iname %in% duplicated$pert_iname &
                                  pert_type %in% duplicated$pert_type &
                                  cell_id %in% duplicated$cell_id &
                                  pert_time %in% duplicated$pert_time)]
    PerturbAnno = rbind(PerturbAnno, duplicated_IDs)
  }

  PerturbAnno
}
