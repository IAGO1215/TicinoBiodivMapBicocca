### This R script is aimed at calculating the alpha and beta biodiversity indices from the Sentinel-2 raster data ###
### The workflow is: import raster -> calculate PCA and select PCs -> calculate spectral species -> calculate alpha and beta ###

# Clean environment
rm(list=ls(all=TRUE));gc()
# Import the R package
library(biodivMapR)
# Set absolute working directory
path_abs <- "C:/Users/m1865/Desktop/Ticino"
setwd(path_abs)

### Reading files ###

print("----LOADING INPUT SATELLITE IMAGE FILES----")

# Importing raster files
# The raw data folder
dir_prodata <- file.path(path_abs, 'ProcessedData')
# Get the names of all the subfolders
name_prodata <- list.dirs(dir_prodata, full.name = FALSE, recursive = FALSE)
# Get all the subfolders in our "ProcessedData" folder
subdir_prodata <- list.dirs(dir_prodata, recursive = FALSE)
# Get the paths to all the raster files! 
path_raster <- list()
for (x in 1:length(subdir_prodata)) {
  temp <- file.path(subdir_prodata[x], paste0(name_prodata[x],'Cropped'))
  print(temp)
  path_raster <- append(path_raster, temp)
}
# Mask paths
path_mask <- list()
for (x in 1:length(subdir_prodata)) {
  temp <- file.path(subdir_prodata[x], paste0(name_prodata[x],'Mask'))
  print(temp)
  path_mask <- append(path_mask, temp)
}
# Output directory (No need to create the subfolders on our own since biodiverR will create them automatically)
outDir <- file.path(path_abs,'Results')

### Setting parameters ###

print("----SETTING PARAMETERS FOR BIODIVERSITY----")

# continuum removal
continuum_Removal <- TRUE
# Type of dimensionality reduction
typePCA <- 'SPCA'
# Automatically set to FALSE if TypePCA     = 'MNF'
filterPCA <- FALSE
# window size forcomputation of spectral diversity
window_size <- 10
# computational parameters
nbCPU <- 4
MaxRAM <- 0.2
# number of clusters (spectral species)
nbClusters <- 50
# spectral filtering
spectral_excluded <- NULL

### PCA ###

print("----PERFORMING PCA----")

list_PCA <- list()
for (x in 1:length(subdir_prodata)){
  PCA_Output <- biodivMapR::perform_PCA(Input_Image_File = path_raster[[x]],
                                        Output_Dir = outDir,
                                        Input_Mask_File = path_mask[[x]],
                                        TypePCA = typePCA,
                                        FilterPCA = filterPCA,
                                        nbCPU = nbCPU,
                                        MaxRAM = MaxRAM,
                                        Continuum_Removal = continuum_Removal, 
                                        Excluded_WL = spectral_excluded)
  print(paste0('PCA generated for ',path_raster[[x]]))
  list_PCA <- append(list_PCA, PCA_Output)
}

print("----PERFORMING PCA SUCCESSFULLY----")

print("----MANUALLY SELECTING PCs----")
selected_PCs <- c(1,2,3,4)