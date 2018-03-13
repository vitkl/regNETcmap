##' Read Connectivity Map data given perturbation annotations
##' @rdname readCMAP
##' @name readCMAP
##' @author Vitalii Kleshchevnikov
##' @description \code{readCMAP} reads Connectivity Map data given perturbation annotations for a set of HUGO gene names (or Entrez gene ID) including additional filtering by cell line, perturbation type (compound, shRNA, overexpression, e.g.) and time. Details \code{\link{geneName2PerturbAnno}}
##' @param PerturbAnno a character vector of HUGO gene names
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
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
                    CMap_files) {
  unzipped = substr(CMap_files$sig_level5[2],
                    1, nchar(CMap_files$sig_level5[2])-3)
  if(!file.exists(unzipped)) unzipCMapData(CMap_files)
  CMAP = parse.gctx(fname = unzipped, rid = NULL, cid = PerturbAnno$sig_id, set_annot_rownames = F,
             matrix_only = F)
  CMAP = annotate.gct(CMAP, PerturbAnno, dim="col", keyfield="sig_id")
  CMAP = annotate.gct(CMAP, openFeatureData(CMap_files), dim="row", keyfield="pr_gene_id")
  CMAP
}
