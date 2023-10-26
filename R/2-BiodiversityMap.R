### This R script is aimed at calculating the alpha and beta biodiversity indices from the Sentinel-2 raster data ###
### The workflow is: import rasters in different subfolders -> calculate PCA and select PCs -> calculate spectral species -> calculate alpha and beta ###

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
  list_PCA <- append(list_PCA, list(PCA_Output))
}

print("----PERFORMING PCA SUCCESSFULLY----")

print("----MANUALLY SELECTING PCs----")
selected_PCs <- c(1,2,3,4)

## TTTTTTTTTTEST
library(raster)
mask1 <- raster(path_mask[[1]])
Mask <- stars::read_stars(path_mask[[1]], proxy = F)
ncols <- as.integer(ncol(mask1))
nrows <- as.integer(nrow(mask1))
HDR <- read_ENVI_header(get_HDR_name(path_raster[[1]]))
!ncols == HDR$samples | !nrows == HDR$lines
ncols == HDR$samples & nrows == HDR$lines

### Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))
list_spectral <- list()
for (x in 1:length(subdir_prodata)){
  Kmeans_info <- biodivMapR::map_spectral_species_py(Input_Image_File = path_raster[[x]],
                                                  Input_Mask_File = path_mask[[x]],
                                                  Output_Dir = outDir,
                                                  SpectralSpace_Output = list_PCA[[x]],
                                                  SelectedPCs = selected_PCs,
                                                  nbclusters = nbClusters,
                                                  nbCPU = nbCPU, 
                                                  MaxRAM = MaxRAM)
  print(paste0('Spectral species generated for ',path_raster[[x]]))
  list_spectral <- append(list_spectral, list(Kmeans_info))
}

print("----SPECTRAL SPECIES COMPUTED----")

### Alpha Beta Diversity ###
print("----MAP ALPHA DIVERSITY----")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')

for (x in 1:length(subdir_prodata)){
  biodivMapR::map_alpha_div(Input_Image_File = path_raster[[x]],
                Input_Mask_File = path_mask[[x]],
                Output_Dir = outDir,
                TypePCA = typePCA,
                window_size = window_size,
                nbCPU = nbCPU,
                MaxRAM = MaxRAM,
                Index_Alpha = Index_Alpha,
                nbclusters = nbClusters)
  print(paste0('Alpha diversity map generated for ',path_raster[[x]]))
}

print("----MAP BETA DIVERSITY----")
for (x in 1:length(subdir_prodata)){
  biodivMapR::map_beta_div(Input_Image_File = path_raster[[x]],
               Output_Dir = outDir,
               TypePCA = typePCA,
               window_size = window_size,
               nbCPU = nbCPU,
               MaxRAM = MaxRAM,
               nbclusters = nbClusters)
  print(paste0('Beta diversity map generated for ',path_raster[[x]]))
}

print("----BETA DIVERSITY COMPUTED----")

### Field Splot Biodiversity ###

print("----FIELD PLOT STEP----")
# location of the directory where shapefiles used for validation are saved
vector_Dir <- file.path(path_abs,"FieldData","FieldPlotsPolygon")
# list vector data (In our case, there is only one shapefile, so alternatively we can directly input its abs path)
vector_Path <- biodivMapR::list_shp(vector_Dir)
vector_Name <- tools::file_path_sans_ext(basename(vector_Path))

Biodiv_Indicators <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral$SpectralSpecies, 
                                          Plots = vector_Path,
                                          nbclusters = nbClusters, 
                                          Raster_Functional = list_PCA[[1]]$PCA_Files, 
                                          Selected_Features = selected_PCs)
