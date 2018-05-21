##' Calculate and remove principal components from Connectivity map data
##' @rdname svdCMAP
##' @name reconstructCMAP
##' @description \code{reconstructCMAP}: calculate and remove principal components from Connectivity map data. If CMAP data matrix have already been reconstructed using the same PCs and pert_types this function will read previously saved SVD output and reconstructed CMAP data matrix
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param remove_PCs integer or character, indexes - which principal components to remove. "none" to keep all PCs
##' @param return_PCs logical, return the output of \link[base]{svd} decomposition
##' @param return_CMAP logical, return Connectivity map data reconstructed using only specified PCs
##' @param save_CMAP logical, save Connectivity map data reconstructed using only specified PCs. Returned \code{CMap_files} contain the path to reconstructed Connectivity map data matrix and corresponding supplementary files.
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}. Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#
##' @return \code{reconstructCMAP}:
##' @import data.table
##' @export reconstructCMAP
##' @examples
##' library(regNETcmap)
##' CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
reconstructCMAP = function(CMap_files, remove_PCs = "none", return_PCs = T, return_CMAP = T, save_CMAP = F, pert_types = c("trt_sh.cgs", "trt_oe", "trt_cp", "trt_lig", "trt_oe.mut", "trt_xpr", "trt_sh.css", "ctl_vehicle.cns", "ctl_vector.cns", "ctl_untrt.cns")){

  types = perturbTable(CMap_files, ~ pert_type)
  pert_types = pert_types[pert_types %in% names(types)]

  filename = paste0("res_",
                    digest::sha1(x = paste0(c(remove_PCs,
                                              return_PCs, return_CMAP, save_CMAP,
                                              pert_types))))
  filepath = R.utils::filePath(dirname(CMap_files$sig[2]), filename)
  filepath.gctx = paste0(filepath, ".gctx")
  filepath.svd = paste0(filepath, ".svd")

  if(file.exists(filepath)){
    load(filepath)
    svd_obj = res$svd_obj
    cmap = res$cmap
    CMap_files = res$CMap_files
  } else {

    cmap = readCMAPsubset(is_touchstone = "all", pert_types = pert_types, pert_itimes = c("all"), cell_ids = c("all"), CMap_files)

    svd_obj = svd(cmap@mat)

    if(remove_PCs[1] != "none") {
      svd_obj$DtV = tcrossprod(diag(svd_obj$d)[-remove_PCs,-remove_PCs], svd_obj$v[,-remove_PCs])
      Yhat <- svd_obj$u[,-remove_PCs] %*% svd_obj$DtV
    } else {
      svd_obj$DtV = tcrossprod(diag(svd_obj$d), svd_obj$v)
      Yhat <- svd_obj$u %*% svd_obj$DtV
    }

    colnames(Yhat) = colnames(cmap@mat)
    rownames(Yhat) = rownames(cmap@mat)
    cmap@mat = Yhat
    rm(Yhat)
    if(save_CMAP){
      cmapR::write.gctx(ds = cmap, ofile = filepath.gctx,
                        appenddim = T, compression_level = 0,
                        matrix_only = F, max_chunk_kb = 1024)
      CMap_files$sig[2] = filepath.gctx
    }
  }

  if(!return_CMAP) cmap = NA
  if(!return_PCs) svd_obj = NA
  res = list(svd_obj = svd_obj,
             cmap = cmap,
             CMap_files = CMap_files)

  save(res, file = filepath)
  res
}

##' Calculate principal components (SVD decomposition) from Connectivity map data
##' @rdname svdCMAP
##' @name svdCMAP
##' @description \code{svdCMAP}: calculate and remove principal components from Connectivity map data. If CMAP data matrix have already been reconstructed using the same PCs and pert_types this function will read previously saved SVD output and reconstructed CMAP data matrix
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param return_PCs logical, return the output of \link[base]{svd} decomposition. Alternatively just calculate and save.
##' @param is_touchstone logical, select only the perturbation data that is a part of the touchstone dataset (genes profiled across all 9 core cell lines)? Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}. Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#
##' @param pert_itimes a character vector of perturbation times, query for available perturbation times using \code{perturbTable(CMap_files, ~ pert_itime)}
##' @param cell_ids a character vector of cell line names, query for available cell line names using \code{perturbTable(CMap_files, ~ cell_id)}
##' @param gene_names a character vector of HUGO gene names. Use \code{"all"} to select all perturbations.
##' @param keep_one_oe keep only one overexpression experiment per gene and condition. Applicable only when \code{pert_types = "trt_oe"}. Perturbations where sig_id matches pert_id are retained ("one"). To invert the selection use "other". To select all use "all"
##' @param landmark_only read only landmark genes (rows, measured). Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit
##' @param ... passed to \link[base]{svd}: \code{nu} and \code{nv}
##' @return \code{svdCMAP}:
##' @import data.table
##' @export svdCMAP
##' @examples
##' library(regNETcmap)
##' CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
svdCMAP = function(CMap_files, return_PCs = T, is_touchstone = "all",
                   pert_types = c("trt_sh.cgs", "trt_oe", "trt_cp", "trt_lig", "trt_oe.mut",
                                  "trt_xpr", "trt_sh.css", "ctl_vehicle.cns", "ctl_vector.cns", "ctl_untrt.cns"),
                   pert_itimes = c("all"), cell_ids = c("all"), gene_names = c("all"),
                   keep_one_oe = "all", landmark_only = F, ...){

  types = perturbTable(CMap_files, ~ pert_type)
  pert_types = pert_types[pert_types %in% names(types)]

  filename = paste0("res_",
                    digest::sha1(x = paste0(c(is_touchstone, pert_types, pert_itimes,
                                              cell_ids, gene_names, keep_one_oe, landmark_only))))
  filepath = R.utils::filePath(dirname(CMap_files$sig[2]), filename)
  filepath.svd = paste0(filepath, ".svd")

  if(file.exists(filepath.svd)){
    load(filepath)
  } else {
    cmap = readCMAPsubset(is_touchstone = is_touchstone, pert_types = pert_types,
                          pert_itimes = pert_itimes, cell_ids = cell_ids, CMap_files,
                          gene_names = gene_names, keep_one_oe = keep_one_oe,
                          landmark_only = landmark_only)

    svd_obj = svd(cmap@mat, ...)

    colnames_cmap = colnames(cmap@mat)
    rownames_cmap = rownames(cmap@mat)
    svd_obj$cdesc = cmap@cdesc
    svd_obj$rdesc = cmap@rdesc
    rm(cmap)

    svd_obj$DtV = tcrossprod(diag(svd_obj$d), svd_obj$v)

    rownames(svd_obj$v) = colnames_cmap
    rownames(svd_obj$u) = rownames_cmap

    save(svd_obj, file = filepath.svd)
    if(return_PCs) return(svd_obj) else message("svd decomposition was already performed on this dataset: use return_PCs = T to obtain results")
  }
}
