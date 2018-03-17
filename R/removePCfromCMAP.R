##' Calculate and remove principal components from Connectivity map data
##' @rdname removePCfromCMAP
##' @name removePCfromCMAP
##' @description \code{removePCfromCMAP}: calculate and remove principal components from Connectivity map data. If CMAP data matrix have already been reconstructed using the same PCs and pert_types this function will read previously saved SVD output and reconstructed CMAP data matrix
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @param keep_PCs integer, indexes - which principal components to keep
##' @param return_PCs logical, return the output of \link[base]{svd} decomposition
##' @param return_CMAP logical, return Connectivity map data reconstructed using only specified PCs
##' @param save_CMAP logical, save Connectivity map data reconstructed using only specified PCs. Returned \code{CMap_files} contain the path to reconstructed Connectivity map data matrix and corresponding supplementary files.
##' @param pert_types a character vector of perturbation types, query for available perturbation types using \code{perturbTable(CMap_files, ~ pert_type)}. Details: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#
##' @return \code{removePCfromCMAP}: data.table containing gene set id (single), which genes may regulate this set, test statistic and difference in medians between each gene set and all other genes
##' @import data.table
##' @export removePCfromCMAP
removePCfromCMAP = function(CMap_files, keep_PCs = "all", return_PCs = T, return_CMAP = T, save_CMAP = F, pert_types = c("trt_sh.cgs", "trt_oe", "trt_cp", "trt_lig", "trt_oe.mut", "trt_xpr", "trt_sh.css", "ctl_vehicle.cns", "ctl_vector.cns", "ctl_untrt.cns")){

  types = perturbTable(CMap_files, ~ pert_type)
  pert_types = pert_types[pert_types %in% names(types)]

  filename = paste0("res_", digest::sha1(x = paste0(c(keep_PCs, pert_types))))
  filepath = R.utils::filePath(dirname(CMap_files$sig_level5[2]), filename)
  filepath.gctx =paste0(filepath, ".gctx")

  if(file.exists(filepath)){
    load(filepath)
    svd_obj = res$svd_obj
    cmap = res$cmap
    CMap_files = res$CMap_files
  } else {

    if(save_CMAP){
      cmapR::write.gctx(ds = cmap, ofile = filepath.gctx,
                        appenddim = T, compression_level = 0,
                        matrix_only = F, max_chunk_kb = 1024)
      CMap_files$sig_level5[2] = filepath.gctx
    }
  }

  res = list(svd_obj = svd_obj,
             cmap = cmap,
             CMap_files = CMap_files)
  if(!return_CMAP) res$cmap = NA
  if(!return_PCs) res$svd_obj = NA
  save(res, file = filepath)
  res
}
