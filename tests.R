library(roxygen2)
library(devtools)
document()
install()

library(regNETcmap)
?loadCMap

CMap_dirs = loadCMap(directory = "./test_data/")
fdata = openFeatureData(CMap_dirs)
pdata = openpData(CMap_dirs)
perturbTable(CMap_dirs)
perturbTable(CMap_dirs, ~ pert_time)

p53IDs = geneName2PertID(CMap_dirs = CMap_dirs)

pertTypesTimes(CMap_dirs, p53IDs)
