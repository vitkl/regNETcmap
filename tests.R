library(roxygen2)
library(devtools)
document()
install()

library(regNETcmap)
?loadCMap

CMap_dirs = loadCMap(directory = "./test_data/")
fdata = openFeatureData(CMap_dirs)
pdata = openpData(CMap_dirs)

pertTypes = function(CMap_dirs){
  table(openPerturbDetails(CMap_dirs)$pert_type)
}

pertTypes(CMap_dirs)

pertTypesTimes = function(CMap_dirs, pert_ids = NULL, pdata = NULL){

}

pertTypesTimes(CMap_dirs)

pertTable = function(CMap_dirs, formula = pert_type ~ pert_time, pert_ids = NULL, pdata = NULL) {
  if(is.null(pdata)) pdata = openpData(CMap_dirs)
  if(is.null(pert_ids)) {
    table(pdata[,.(eval(formula[[2]]), eval(formula[[3]]))],
          dnn = c(as.character(formula[[2]]), as.character(formula[[3]])))
  } else {
    table(pdata[distil_id %in% pert_ids$distil_id,.(eval(formula[[2]]), eval(formula[[3]]))],
          dnn = c(as.character(formula[[2]]), as.character(formula[[3]])))
  }
}
pertTable(CMap_dirs)

geneName2PertID = function(gene_names = "TP53", CMap_dirs, is_touchstone = c("all", T, F), pert_types = c("trt_sh.cgs", "trt_sh", "trt_sh.css"), pert_times) {
  pdata = openpData(CMap_dirs)[pert_type %in% pert_types & pert_iname %in% gene_names]
  pdetails = openPerturbDetails(CMap_dirs)[,.(pert_id, pert_type, is_touchstone)][pert_type %in% pert_types]
  pIDdata = pdata[pdetails, on = c("pert_id", "pert_type"), nomatch = 0]

  if(is_touchstone == "all") pIDdata = pIDdata else {
    if(is_touchstone) pIDdata = pIDdata[is_touchstone == 1] else if(!is_touchstone) pIDdata = pIDdata[is_touchstone == 0]
  }
  pIDdata
}

p53IDs = geneName2PertID(CMap_dirs = CMap_dirs)

pertTypesTimes(CMap_dirs, p53IDs)
