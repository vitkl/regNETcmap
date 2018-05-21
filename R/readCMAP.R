##' Read Connectivity Map data given perturbation annotations
##' @rdname readCMAP
##' @name readCMAP
##' @author Vitalii Kleshchevnikov
##' @description \code{readCMAP} reads Connectivity Map data given perturbation annotations for a set of HUGO gene names (or Entrez gene ID) including additional filtering by cell line, perturbation type (compound, shRNA, overexpression, e.g.) and time. Details \code{\link{geneName2PerturbAnno}}
##' @param PerturbAnno data.table containing the perturbation details from the Connectivity map project produced by \code{\link{geneName2PerturbAnno}}
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param landmark_only read only landmark genes (rows, measured). Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit
##' @return object of class 'GCT' [package "cmapR"] with 7 slots containing z-score matrix, perturbation and feature details of the Connectivity map project
##' @export readCMAP
##' @seealso \code{\link{openCellInfo}}, \code{\link{loadCMap}}, \code{\link{perturbTable}}
##' @examples
##' # read P53 knockdown perturbations
##' CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
##' p53IDs = geneName2PerturbAnno(CMap_files = CMap_files, is_touchstone = T, pert_types = "trt_sh.cgs")
##' CMAP = readCMAP(PerturbAnno = p53IDs, CMap_files = CMap_files)
##' # read all knockdown perturbations (3.4GB)
##' allIDs = geneName2PerturbAnno(gene_names = "all", CMap_files = CMap_files, is_touchstone = T, pert_types = "trt_sh.cgs")
##' CMAP = readCMAP(PerturbAnno = allIDs, CMap_files = CMap_files)
readCMAP = function(PerturbAnno,
                    CMap_files,
                    landmark_only = F) {
  unzipped = substr(CMap_files$sig[2],
                    1, nchar(CMap_files$sig[2])-3)
  if(!file.exists(unzipped)) unzipCMapData(CMap_files)
  fdata = openFeatureData(CMap_files)
  if(landmark_only) fdata = fdata[pr_is_lm == 1]
  if(CMap_files$level == 4) key_field = "inst_id"
  if(CMap_files$level == 5) key_field = "sig_id"
  CMAP = parse.gctx(fname = unzipped, rid = fdata$pr_gene_id, cid = unlist(PerturbAnno[,key_field,with = F]), set_annot_rownames = F,
             matrix_only = F)
  if(ncol(CMAP@cdesc) == 1) CMAP = annotate.gct(CMAP, PerturbAnno, dim="col", keyfield=key_field)
  if(ncol(CMAP@rdesc) == 1) CMAP = annotate.gct(CMAP, fdata, dim="row", keyfield="pr_gene_id")
  CMAP
}
