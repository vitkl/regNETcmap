#devtools::install_github("klutometis/roxygen")
#devtools::install_github("cmap/cmapR")
library(roxygen2)
library(devtools)
document()
install()

devtools::install_github("vitkl/regNETcmap", dependencies = T)

#devtools::install_github("cmap/cmapR")
#devtools::install_github("vitkl/regNETcmap")
#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library(regNETcmap)
?loadCMap

CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
fdata = openFeatureData(CMap_files)
pdata = openpData(CMap_files)
perturbTable(CMap_files)
perturbTable(CMap_files, pert_type ~ cell_id)
perturbTable(CMap_files, ~ pert_time)

p53IDs = geneName2PerturbAnno(gene_names = "TP53", CMap_files,
                              is_touchstone = c("all", T, F), pert_types = c("trt_sh.cgs"), pert_times = c("all"), cell_ids = c("all"))

pertTypesTimes(CMap_files, p53IDs)

pdata[pert_type %in% "trt_sh", uniqueN(pert_iname)]
pdata[pert_type %in% "trt_oe", uniqueN(pert_iname)]
pdata[pert_type %in% "trt_sh.cgs", uniqueN(pert_iname)]
pdata[pert_type %in% c("trt_sh", "trt_sh.cgs"), uniqueN(pert_iname)]

uniqueN(intersect(pdata[pert_type %in% "trt_oe", pert_iname],
                  pdata[pert_type %in% "trt_sh.cgs", pert_iname]))

pdata[pert_type %in% "trt_cp", uniqueN(pert_iname)]
pdata[pert_type %in% "trt_sh.css", uniqueN(pert_iname)]

par(mar = c(5,5,5,0), bg = "black", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white")
barplot(c(50, 10, 12, 40), col = c("#EC2E53", "#328BFF","#EC2E53", "#328BFF"),
        names.arg = c("same phenotype\nactivating", "same phenotype\ninhibiting","other phenotype\nactivating", "other phenotype\ninhibiting"),
        main = "number of genes regulating transcriptional phenotype",
        cex.axis = 2, cex.main = 2, cex.names = 1.71, border = "transparent", mgp =  c(3, 2.3, -0.2))

ks.test_expl(alternative = "less", prop_pos = 0.05, prop_neg = 0.25, seed = 1)

# read all knockdown perturbations (3.4GB)
CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
allIDs = geneName2PerturbAnno(gene_names = "all", CMap_files = CMap_files,
                              is_touchstone = T, pert_types = "trt_oe",
                              pert_times = 96, cell_ids = "A375")
perturbTable(pdata = allIDs)
duplicated = table(allIDs$pert_iname)[table(allIDs$pert_iname) > 1]
duplicated
duplicated_IDs = geneName2PerturbAnno(gene_names = names(duplicated), CMap_files = CMap_files,
                                      is_touchstone = T, pert_types = "trt_oe",
                                      pert_times = 96, cell_ids = "A375")
duplicated_IDs = duplicated_IDs[order(pert_iname)]
duplicated_IDs[, pert_id_in_sig_id := grep(pert_id, sig_id, ignore.case = T, value = T), by = pert_id]
duplicated_IDs = duplicated_IDs[sig_id == pert_id_in_sig_id]

CMAP = readCMAP(PerturbAnno = allIDs, CMap_files = CMap_files)

CMAPsub = readCMAPsubset(is_touchstone = T, pert_types = "trt_sh.cgs",
                         pert_times = 96, cell_ids = "A375", CMap_files)
perturbTable(pdata = as.data.table(CMAPsub@cdesc))

norm_counts = assays(data)$norm_counts
norm_counts[norm_counts < 1] = NA
row_var = rowVars(norm_counts, na.rm = T)
hist(log10(norm_counts[order(row_var, decreasing = T)[1],]), main = 1)

for (i in 1:100) {
  hist(log10(norm_counts[order(row_var, decreasing = T)[i],]), main = i)
  Sys.sleep(1)
}


# cluster means
norm_counts = assays(data)$norm_counts
norm_counts_clusters = unique(clusters)
row_med = sapply(norm_counts_clusters, function(cluster, norm_counts) {
  rowMedians(norm_counts[,grepl(cluster,colnames(norm_counts))])
}, norm_counts)
row_med = row_med[rowSums2(row_med != 0) > 1,]
row_var = rowVars(row_med, na.rm = T)
hist(log10(row_med[order(row_var, decreasing = T)[1],]), main = 1)
plot(density(log10(row_med[order(row_var, decreasing = T)[1],])), main = 1)

for (i in 1:100) {
  #hist(log10(row_med[order(row_var, decreasing = T)[i],]), main = i)
  plot(density(log10(row_med[order(row_var, decreasing = T)[i],])), main = i)
  Sys.sleep(1)
}

STAT3_OE_PC3 = readCMAPsubset(is_touchstone = T,
                              pert_types = "trt_oe",
               pert_times = 96,
               cell_ids = "PC3", CMap_files,
               gene_names = "STAT3", keep_one_oe = "one",
               landmark_only = F)
STAT3_OE_PC3@mat[rownames(STAT3_OE_PC3@mat) %in% "6774",]
STAT3_OE_PC3@mat[STAT3_OE_PC3@rdesc$pr_gene_symbol %in% "STAT3",]

STAT1_OE_PC3 = readCMAPsubset(is_touchstone = T,
                              pert_types = "trt_oe",
                              pert_times = 96,
                              cell_ids = "PC3", CMap_files,
                              gene_names = "STAT1", keep_one_oe = "one",
                              landmark_only = F)
STAT1_OE_PC3@mat[rownames(STAT1_OE_PC3@mat) %in% "6772",]
STAT1_OE_PC3@mat[STAT1_OE_PC3@rdesc$pr_gene_symbol %in% "STAT1",]

MYC_SH_HT29 = readCMAPsubset(is_touchstone = T,
                              pert_types = "trt_sh.cgs",
                              pert_times = 96,
                              cell_ids = "HT29", CMap_files,
                              gene_names = "MYC", keep_one_oe = "one",
                              landmark_only = F)
MYC_SH_HT29@mat[rownames(MYC_SH_HT29@mat) %in% "4609",]
MYC_SH_HT29@mat[MYC_SH_HT29@rdesc$pr_gene_symbol %in% "MYC",]

STAT1_SH_HT29 = readCMAPsubset(is_touchstone = T,
                             pert_types = "trt_sh.cgs",
                             pert_times = 96,
                             cell_ids = "HT29", CMap_files,
                             gene_names = "STAT1", keep_one_oe = "one",
                             landmark_only = F)
STAT1_SH_HT29@mat[rownames(STAT1_SH_HT29@mat) %in% "6772",]
STAT1_SH_HT29@mat[STAT1_SH_HT29@rdesc$pr_gene_symbol %in% "STAT1",]

control_SH_HT29 = readCMAPsubset(is_touchstone = F,
                               pert_types = "ctl_vector.cns",
                               pert_times = 96,
                               cell_ids = "HT29", CMap_files,
                               gene_names = "EMPTY_VECTOR", keep_one_oe = "one",
                               landmark_only = F)
control_SH_HT29@mat[rownames(control_SH_HT29@mat) %in% "6772",]
control_SH_HT29@mat[control_SH_HT29@rdesc$pr_gene_symbol %in% "STAT1",]

## phase 1 vs phase 2##
# phase 1
CMap_files = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/",
                      level = 4, phase = 1)
fdata = openFeatureData(CMap_files)
pdata = openpData(CMap_files)
cell_info = openCellInfo(CMap_files)
perturb_details = openPerturbDetails(CMap_files)

# phase 2
CMap_files2 = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/",
                      level = 4, phase = 2)
fdata2 = openFeatureData(CMap_files2)
cell_info2 = openCellInfo(CMap_files2)
pdata2 = openpData(CMap_files2)
perturb_details2 = openPerturbDetails(CMap_files2)
