##' Get Connectivity Map perturbation annotations for genes of interest, cell line and perturbation types and times
##' @rdname geneName2PerturbAnno
##' @name geneName2PerturbAnno
##' @author Vitalii Kleshchevnikov
##' @description \code{geneName2PerturbAnno} retrieves perturbation annotations for a set of HUGO gene names (or Entrez gene ID) including additional filtering by cell line, perturbation type (compound, shRNA, overexpression, e.g.) and time
##' @details is_touchstone: A boolean indicating whether the corresponding signature or perturbagen is a member of the Touchstone dataset. Touchstone is a term applied to the subset of CMap perturbagens that are well-annotated and that were systematically profiled across the majority of the core set of 9 cell lines at standardized conditions. Because of these properties Touchstone dataset well-suited as a reference compendium against with to compare external queries.
##' @param gene_names a character vector of HUGO gene names. Use \code{"all"} to select all perturbations.
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param is_touchstone logical, select only the perturbation data that is a part of the touchstone dataset (genes profiled across all 9 core cell lines)? Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}
##' @param pert_times a character vector of perturbation times, query for available perturbation times using \code{perturbTable(CMap_files, ~ pert_time)}
##' @param cell_ids a character vector of cell line names, query for available cell line names using \code{perturbTable(CMap_files, ~ cell_id)}
##' @return data.table containing the perturbation details from the Connectivity map project
##' @import data.table
##' @export geneName2PerturbAnno
##' @seealso \code{\link{openCellInfo}}, \code{\link{loadCMap}}, \code{\link{perturbTable}}
geneName2PerturbAnno = function(gene_names = "TP53", CMap_files, is_touchstone = c("all", T, F), pert_types = c("trt_sh.cgs", "trt_oe"), pert_times = c("all"), cell_ids = c("all")) {
  pdata = openpData(CMap_files)[pert_type %in% pert_types]
  if(!gene_names[1] == "all") pdata = pdata[pert_iname %in% gene_names]
  pdetails = openPerturbDetails(CMap_files)[,.(pert_id, pert_type, is_touchstone)][pert_type %in% pert_types]
  PerturbAnno = pdata[pdetails, on = c("pert_id", "pert_type"), nomatch = 0]

  if(is_touchstone == "all") NULL else {
    if(is_touchstone) PerturbAnno = PerturbAnno[is_touchstone == 1] else if(!is_touchstone) PerturbAnno = PerturbAnno[is_touchstone == 0]
  }
  if(cell_ids[1] == "all") NULL else {
    PerturbAnno = PerturbAnno[cell_id %in% cell_ids]
  }
  if(pert_times[1] == "all") NULL else {
    PerturbAnno = PerturbAnno[pert_time %in% pert_times]
  }

  PerturbAnno
}

##' @rdname geneName2PerturbAnno
##' @name anyGeneID2PerturbAnno
##' @author Vitalii Kleshchevnikov
##' @description \code{anyGeneID2PerturbAnno} retrieves perturbation annotations for a set of gene idenfified using any type of ID in Homo.sapiens or UniProt.ws, including additional filtering by cell line, perturbation type (compound, shRNA, overexpression, e.g.) and time
##' @param keys a character vector of identifiers, check possible options for type of identifier by calling keytypes(org.db)
##' @param org.db org.db of choice: Homo.sapiens or UniProt.ws::UniProt.ws(9606)
##' @param keytype type of identifier used in \code{keys}
##' @return data.table containing the perturbation details from the Connectivity map project
##' @import data.table
##' @export anyGeneID2PerturbAnno
##' @seealso \code{\link{openCellInfo}}, \code{\link{loadCMap}}, \code{\link{perturbTable}}
anyGeneID2PerturbAnno = function(keys = "Q9NZL3", CMap_files, org.db = c("UniProt.ws", "Homo.sapiens")[1], keytype = "UNIPROTKB", is_touchstone = c("all", T, F), pert_types = c("trt_sh.cgs", "trt_sh", "trt_sh.css"), pert_times = c("all"), cell_ids = c("all")){
  if(org.db == "UniProt.ws"){
    db = UniProt.ws::UniProt.ws(taxId = 9606)
    mapping = UniProt.ws::select(x = db, keys = keys, columns = c("GENES", keytype[1]), keytype = keytype[1])
    mapping$GENES = gsub(" .+$", "", mapping$GENES)
    mapping$pert_iname = mapping$GENES
    mapping$GENES = NULL
  }
  if(org.db == "Homo.sapiens"){
    db = Homo.sapiens::Homo.sapiens
    mapping = AnnotationDbi::select(x = db, keys = keys, columns = c("SYMBOL", keytype[1]), keytype = keytype[1])
    mapping$pert_iname = mapping$SYMBOL
    mapping$SYMBOL = NULL
  }

  res = geneName2PerturbAnno(gene_names = unique(mapping$pert_iname), CMap_files = CMap_files, is_touchstone = is_touchstone, pert_types = pert_types, pert_times = pert_times, cell_ids = cell_ids)
  res = res[mapping, on = "pert_iname", nomatch = 0]
  res
}
