##' @name writeLandMarkOnlyFile
##' @rdname loadCMap
##' @author Vitalii Kleshchevnikov
##' @description \code{writeLandMarkOnlyFile()} generates and writes .gctx file containing only landmark genes from corresponding Connectivity Map release. This requires about 32 GB of RAM.
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @export writeLandMarkOnlyFile
writeLandMarkOnlyFile = function(CMap_files = loadCMap(directory = "./data/cmap/",
                                                    level = 4, phase = 1)){
  landmark_only = T
  pdata = openpData(CMap_files)
  pdetails = openPerturbDetails(CMap_files)[,.(pert_id, pert_type, is_touchstone)]
  PerturbAnno = merge(pdata, pdetails, by = c("pert_id", "pert_type"), all.x = T, all.y = F)
  PerturbAnno[, is_touchstone := as.numeric(is_touchstone)]
  PerturbAnno[is.na(is_touchstone), is_touchstone := 0]
  PerturbAnno[, is_touchstone := as.factor(is_touchstone)]

  fdata = openFeatureData(CMap_files)
  unzipped = substr(CMap_files$sig[2],
                    1, nchar(CMap_files$sig[2])-3)
  if(!file.exists(unzipped)) unzipCMapData(CMap_files)
  if(landmark_only) fdata = fdata[pr_is_lm == 1]
  if(CMap_files$level == 4) key_field = "inst_id"
  if(CMap_files$level == 5) key_field = "sig_id"
  CMAP = parse.gctx(fname = unzipped, rid = fdata$pr_gene_id, cid = unlist(PerturbAnno[,key_field,with = F]), set_annot_rownames = F,
                    matrix_only = F)
  if(ncol(CMAP@cdesc) == 1) CMAP = annotate.gct(CMAP, PerturbAnno, dim="col", keyfield=key_field)
  if(ncol(CMAP@rdesc) == 1) CMAP = annotate.gct(CMAP, fdata, dim="row", keyfield="pr_gene_id")

  unzipped = substr(CMap_files$sig[3],
                    1, nchar(CMap_files$sig[3])-3)
  write.gctx(ds = CMAP, ofile = unzipped,
             appenddim = F, compression_level = 0,
             matrix_only = F, max_chunk_kb = 512)
  unzipped
}
