#devtools::install_github("klutometis/roxygen")
#devtools::install_github("cmap/cmapR")
library(roxygen2)
library(devtools)
document()
install()

#devtools::install_github("cmap/cmapR")
#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library(regNETcmap)
?loadCMap

CMap_dirs = loadCMap(directory = "../regulatory_networks_by_cmap/data/cmap/")
fdata = openFeatureData(CMap_dirs)
pdata = openpData(CMap_dirs)
perturbTable(CMap_dirs)
perturbTable(CMap_dirs, pert_type ~ cell_id)
perturbTable(CMap_dirs, ~ pert_time)

p53IDs = geneName2PerturbAnno(gene_names = "TP53", CMap_dirs,
                              is_touchstone = c("all", T, F), pert_types = c("trt_sh.cgs"), pert_times = c("all"), cell_ids = c("all"))

pertTypesTimes(CMap_dirs, p53IDs)

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



