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

CMap_dirs = loadCMap(directory = "./test_data/")
fdata = openFeatureData(CMap_dirs)
pdata = openpData(CMap_dirs)
perturbTable(CMap_dirs)
perturbTable(CMap_dirs, pert_type ~ cell_id)
perturbTable(CMap_dirs, ~ pert_time)

p53IDs = geneName2PertID(CMap_dirs = CMap_dirs)

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
        names.arg = c("same module\nactivating", "same module\ninhibiting","other module\nactivating", "other module\ninhibiting"),
        main = "number of gene-module interactions",
        cex.axis = 2, cex.main = 2, cex.names = 1.71, border = "transparent", mgp =  c(3, 2.3, -0.2))

ks.test_expl(alternative = "greater", prop_pos = 0.33, prop_neg = 0.4, seed = 1)
