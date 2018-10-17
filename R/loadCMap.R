##' Load Connectivity Map project data (1st and 2nd phase, Level 4 and 5)
##' @name loadCMap
##' @author Vitalii Kleshchevnikov
##' @description Download Connectivity Map project data to a specified directory. !! ATTENTION !! This is a 5-20 GB download. Details: \link{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138}, \link{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742}. A description of the latest file, and a table listing the contents of the 'Broad_LINCS_auxiliary_datasets.tar.gz' file are updated in the following document: \link{https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#heading=h.l6bq0r1aih50}
##' @details LINCS aims to enable a functional understanding of biology by cataloging changes in gene expression and other cellular processes that occur when cells are exposed to a variety of perturbing agents. The Broad Institute LINCS Center for Transcriptomics contributes to this collaborative effort by application of the Connectivity Map concept. In brief, the study design involves the generation of a compendium of transcriptional expression data from cultured human cells treated with small-molecule and genetic loss/gain of function perturbagens. These measurements are made using the L1000 high-throughput gene-expression assay that enables data generation at an unprecedented scale. The data are processed through a computational system, that converts raw fluorescence intensities into differential gene expression signatures. The data at each stage of the pre-processing are available:
##' @details Level 1 (LXB) - raw, unprocessed flow cytometry data from Luminex scanners. One LXB file is generated for each well of a 384-well plate, and each file contains a fluorescence intensity value for every observed analyte in the well.
##' @details Level 2 (GEX) - gene expression values per 1,000 genes after deconvolution from Luminex beads.
##' @details Level 3 (Q2NORM) - gene expression profiles of both directly measured landmark transcripts plus inferred genes. Normalized using invariant set scaling followed by quantile normalization.
##' @details Level 4 (Z-SCORES) - signatures with differentially expressed genes computed by robust z-scores for each profile relative to control (PC relative to plate population as control; VC relative to vehicle control).
##' @details Level 5 (SIG) consists of the replicates, usually 3 per treatment, aggregated into a single differential expression vector derived from the weighted averages of the individual replicates.
##' @param directory dir where to save Connectivity Map project data
##' @return list of paths to Connectivity Map project files
##' @import cmapR
##' @import downloader
##' @export loadCMap
loadCMap = function(directory = getwd(), level = 5, phase = 1, landmark_only = T){

  ## phase 1 ##
  if(phase == 1){
    base_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl"
    # level 4
    if(level == 4){
      sig_url = paste0(base_url,"/GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx.gz")
      sig_file = paste0(directory,"/GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx.gz")
      pdata_url = paste0(base_url,"/GSE92742_Broad_LINCS_inst_info.txt.gz")
      pdata_file = paste0(directory,"/GSE92742_Broad_LINCS_inst_info.txt.gz")
    }
    # level 5
    if(level == 5){
      sig_url = paste0(base_url,"/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz")
      sig_file = paste0(directory,"/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz")
      pdata_url = paste0(base_url,"/GSE92742_Broad_LINCS_sig_info.txt.gz")
      pdata_file = paste0(directory,"/GSE92742_Broad_LINCS_sig_info.txt.gz")
    }
    cell_info_url = paste0(base_url,"/GSE92742_Broad_LINCS_cell_info.txt.gz")
    cell_info_file = paste0(directory,"/GSE92742_Broad_LINCS_cell_info.txt.gz")
    featureData_url = paste0(base_url,"/GSE92742_Broad_LINCS_gene_info.txt.gz")
    featureData_file = paste0(directory,"/GSE92742_Broad_LINCS_gene_info.txt.gz")
    perturbation_details_url = paste0(base_url,"/GSE92742_Broad_LINCS_pert_info.txt.gz")
    perturbation_details_file =  paste0(directory,"/GSE92742_Broad_LINCS_pert_info.txt.gz")
  }

  ## phase 2 ##
  if(phase == 2){
    base_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl"
    # level 4
    if(level == 4){
      sig_url = paste0(base_url,"/GSE70138_Broad_LINCS_Level4_ZSPCINF_mlr12k_n345976x12328_2017-03-06.gctx.gz")
      sig_file = paste0(directory,"/GSE70138_Broad_LINCS_Level4_ZSPCINF_mlr12k_n345976x12328_2017-03-06.gctx.gz")
      pdata_url = paste0(base_url,"/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz")
      pdata_file = paste0(directory,"/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz")
    }
    # level 5
    if(level == 5){
      sig_url = paste0(base_url,"/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz")
      sig_file = paste0(directory,"/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz")
      pdata_url = paste0(base_url,"/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
      pdata_file = paste0(directory,"/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
    }
    cell_info_url = paste0(base_url,"/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz")
    cell_info_file = paste0(directory,"/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz")
    featureData_url = paste0(base_url,"/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    featureData_file = paste0(directory,"/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    perturbation_details_url = paste0(base_url,"/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz")
    perturbation_details_file =  paste0(directory,"/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz")
  }

  # generate the name for landmark_only file
  sig_file_landmark_only = basename(sig_file)
  sig_file_landmark_only = gsub("x[[:digit:]]+","x978", sig_file_landmark_only)
  sig_file_landmark_only = gsub("\\.gctx\\.gz","_landmark.gctx.gz", sig_file_landmark_only)
  sig_file_landmark_only = paste0(directory, sig_file_landmark_only)

  CMap_files = list(
    sig = c(sig_url, sig_file, sig_file_landmark_only),
    cell_info = c(cell_info_url, cell_info_file),
    featureData = c(featureData_url, featureData_file),
    pdata = c(pdata_url, pdata_file),
    perturbation_details = c(perturbation_details_url, perturbation_details_file),
    level = level, phase = phase
  )
  if(!dir.exists(directory)) dir.create(directory, recursive = T)
  for(i in 1:(length(CMap_files)-2)){
    # download files only if both unzipped and gzipped versions are not availlable
    zipped = CMap_files[[i]][2]
    unzipped = substr(CMap_files[[i]][2],
                      1, nchar(CMap_files[[i]][2])-3)
    file_not_found = !file.exists(zipped) & !file.exists(unzipped)
    if(i == 1 & landmark_only){
      # do not download data when landmark_only is TRUE and landmark_only data is present, warn and suggest to produce this data when data not present
      zipped = CMap_files[[i]][3]
      unzipped = substr(CMap_files[[i]][3],
                        1, nchar(CMap_files[[i]][3])-3)
      file_not_found = !file.exists(zipped) & !file.exists(unzipped)
      if(file_not_found) {
        warning("Landmark_only file is not found, looking for / downloading full sig file. Following download, use writeLandMarkOnlyFile(CMap_files) to generate and write .gctx file containing only landmark genes. This requires about 32 GB of RAM.")
      } else {
        # change CMap_files[[i]][2] to landmark file so that other functions read the correct file
        CMap_files[[i]][2] = CMap_files[[i]][3]
      }
    }
    if(file_not_found) {
      message(paste0("downloading ", names(CMap_files[i]), " from ", CMap_files[[i]][1]))
      download.file(CMap_files[[i]][1], CMap_files[[i]][2])
    }
  }

  return(CMap_files)
}

##' Read connectivity map annotations
##' @rdname openCellInfo
##' @name openCellInfo
##' @description \code{openCellInfo} reads Connectivity Map project cell line data
##' @author Vitalii Kleshchevnikov
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return data.table containing Cell line Information from the Connectivity map project
##' @importFrom R.utils gunzip
##' @import data.table
##' @export openCellInfo
openCellInfo = function(CMap_files) {
  unzipped = substr(CMap_files$cell_info[2],
                    1, nchar(CMap_files$cell_info[2])-3)
  R.utils::gunzip(CMap_files$cell_info[2],
                  destname = unzipped,
                  remove = F, overwrite = T)
  CellInfo = fread(unzipped, stringsAsFactors = T)
  unlink(unzipped)
  return(CellInfo)
}

##' @rdname openCellInfo
##' @name openFeatureData
##' @author Vitalii Kleshchevnikov
##' @description \code{openFeatureData} reads Connectivity Map project featureData
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return data.table containing featureData (\code{\link[Biobase]{featureData}}) from the Connectivity map project
##' @importFrom R.utils gunzip
##' @import data.table
##' @export openFeatureData
openFeatureData = function(CMap_files) {
  unzipped = substr(CMap_files$featureData[2],
                    1, nchar(CMap_files$featureData[2])-3)
  R.utils::gunzip(CMap_files$featureData[2],
                  destname = unzipped,
                  remove = F, overwrite = T)
  featureData = fread(unzipped, stringsAsFactors = T, colClasses = c("character", "character", "character", "integer", "integer"))
  unlink(unzipped)
  return(featureData)
}

##' @rdname openCellInfo
##' @name openpData
##' @author Vitalii Kleshchevnikov
##' @description \code{openpData} read Connectivity Map project sample descriptions (pData)
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return data.table containing the sample descriptions (pData) (\code{\link[Biobase]{pData}}) from the Connectivity map project
##' @importFrom R.utils gunzip
##' @import data.table
##' @export openpData
openpData = function(CMap_files) {
  unzipped = substr(CMap_files$pdata[2],
                    1, nchar(CMap_files$pdata[2])-3)
  R.utils::gunzip(CMap_files$pdata[2],
                  destname = unzipped,
                  remove = F, overwrite = T)
  pData = fread(unzipped, stringsAsFactors = T, colClasses = c("character"))
  unlink(unzipped)
  return(pData)
}

##' @rdname openCellInfo
##' @name openPerturbDetails
##' @author Vitalii Kleshchevnikov
##' @description \code{openPerturbDetails} read Connectivity Map project perturbation details
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return data.table containing the perturbation details from the Connectivity map project
##' @importFrom R.utils gunzip
##' @import data.table
##' @export openPerturbDetails
openPerturbDetails = function(CMap_files) {
  unzipped = substr(CMap_files$perturbation_details[2],
                    1, nchar(CMap_files$perturbation_details[2])-3)
  R.utils::gunzip(CMap_files$perturbation_details[2],
                  destname = unzipped,
                  remove = F, overwrite = T)

  colClasses = c("character")
  if(CMap_files$level == 5 & CMap_files$phase == 1) colClasses = c("character", "character", "character", "integer", "character", "character", "character", "character")

  PerturbDetails = fread(unzipped, stringsAsFactors = T, colClasses = colClasses)
  unlink(unzipped)
  return(PerturbDetails)
}

##' @rdname loadCMap
##' @name unzipCMapData
##' @author Vitalii Kleshchevnikov
##' @description Unzip Connectivity Map data
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return TRUE
##' @importFrom R.utils gunzip
##' @export unzipCMapData
unzipCMapData = function(CMap_files) {
  unzipped = substr(CMap_files$sig[2],
                    1, nchar(CMap_files$sig[2])-3)
  R.utils::gunzip(CMap_files$sig[2],
                  destname = unzipped,
                  remove = F, overwrite = F)
}

##' @rdname loadCMap
##' @name unzipCMapData
##' @author Vitalii Kleshchevnikov
##' @description Unzip Connectivity Map data
##' @param CMap_files a list of directories and urls produced by \code{\link{loadCMap}}
##' @return TRUE
##' @importFrom R.utils gunzip
##' @export unzipCMapData
