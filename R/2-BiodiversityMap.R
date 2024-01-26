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
dir_Sentinel <- file.path(path_abs, 'ProcessedData')
# Get the names of all the subfolders
name_Sentinel <- list.dirs(dir_Sentinel, full.name = FALSE, recursive = FALSE)
# Get all the subfolders in our "ProcessedData" folder
subdir_Sentinel <- list.dirs(dir_Sentinel, recursive = FALSE)
# Get the paths to all the raster files! 
path_Sentinel_Raster <- list()
for (x in 1:length(subdir_Sentinel)) {
  temp <- file.path(subdir_Sentinel[x], paste0(name_Sentinel[x],'Cropped'))
  print(temp)
  path_Sentinel_Raster <- append(path_Sentinel_Raster, temp)
  rm(temp)
}
# Mask paths
path_Sentinel_Mask <- list()
for (x in 1:length(subdir_Sentinel)) {
  temp <- file.path(subdir_Sentinel[x], paste0(name_Sentinel[x],'Mask'))
  print(temp)
  path_Sentinel_Mask <- append(path_Sentinel_Mask, temp)
  rm(temp)
}
# Output directory (No need to create the subfolders on our own since biodiverR will create them automatically)
outDir_Sentinel <- file.path(path_abs,'Results')

### Setting parameters ###

print("----SETTING PARAMETERS FOR BIODIVERSITY----")

# continuum removal
continuum_Removal <- TRUE
# Type of dimensionality reduction
typePCA <- 'SPCA'
# Automatically set to FALSE if TypePCA     = 'MNF'
filterPCA <- FALSE
# window size forcomputation of spectral diversity
window_size_Sentinel <- 5
# computational parameters
nbCPU <- 4
MaxRAM <- 0.2
# number of clusters (spectral species)
nbClusters_Sentinel <- 20
# spectral filtering - for our case, we don't need any spectrals to be excluded
spectral_excluded <- NULL

### PCA ###

print("----PERFORMING PCA----")

list_PCASentinel <- list()
for (x in 1:length(subdir_Sentinel)){
  temp_PCA <- biodivMapR::perform_PCA(Input_Image_File = path_Sentinel_Raster[[x]],
                                        Output_Dir = outDir_Sentinel,
                                        Input_Mask_File = path_Sentinel_Mask[[x]],
                                        TypePCA = typePCA,
                                        FilterPCA = filterPCA,
                                        nbCPU = nbCPU,
                                        MaxRAM = MaxRAM,
                                        Continuum_Removal = continuum_Removal, 
                                        Excluded_WL = spectral_excluded)
  print(paste0('PCA generated for ',path_Sentinel_Raster[[x]]))
  list_PCASentinel <- append(list_PCASentinel, list(temp_PCA))
  rm(temp_PCA)
}

print("----PERFORMING PCA SUCCESSFULLY----")

print("----MANUALLY SELECTING PCs----")
pc_sel_Sentinel <- c(1,2,3,4,5)

### Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))

list_spectral_Sentinel <- list()
for (x in 1:length(subdir_Sentinel)){
  Kmeans_info <- map_spectral_species_py(Input_Image_File = path_Sentinel_Raster[[x]],
                                                     Input_Mask_File = path_Sentinel_Mask[[x]],
                                                     Output_Dir = outDir_Sentinel,
                                                     SpectralSpace_Output = list_PCASentinel[[x]],
                                                     SelectedPCs = pc_sel_Sentinel,
                                                     nbclusters = nbClusters_Sentinel,
                                                     nbCPU = nbCPU, 
                                                     MaxRAM = MaxRAM)
  print(paste0('Spectral species generated for ',path_Sentinel_Raster[[x]]))
  list_spectral_Sentinel <- append(list_spectral_Sentinel, list(Kmeans_info))
  rm(Kmeans_info)
}

print("----SPECTRAL SPECIES COMPUTED----")

### Alpha Beta Diversity ###
print("----MAP ALPHA DIVERSITY----")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')

for (x in 1:length(subdir_Sentinel)){
  biodivMapR::map_alpha_div(Input_Image_File = path_Sentinel_Raster[[x]],
                Input_Mask_File = path_Sentinel_Mask[[x]],
                Output_Dir = outDir_Sentinel,
                TypePCA = typePCA,
                window_size = window_size_Sentinel,
                nbCPU = nbCPU,
                MaxRAM = MaxRAM,
                Index_Alpha = Index_Alpha,
                nbclusters = nbClusters_Sentinel)
  print(paste0('Alpha diversity map generated for ',path_Sentinel_Raster[[x]]))
}

print("----MAP BETA DIVERSITY----")
for (x in 1:length(subdir_Sentinel)){
  biodivMapR::map_beta_div(Input_Image_File = path_Sentinel_Raster[[x]],
               Output_Dir = outDir_Sentinel,
               TypePCA = typePCA,
               window_size = window_size_Sentinel,
               nbCPU = nbCPU,
               MaxRAM = MaxRAM,
               nbclusters = nbClusters_Sentinel)
  print(paste0('Beta diversity map generated for ',path_Sentinel_Raster[[x]]))
}

print("----BETA DIVERSITY COMPUTED----")

#### Field Plot Biodiversity ####

print("----READ FIELD PLOT SHP----")
# Get the name for all the plot
csv_Plot <- read.csv(file.path(path_abs,"FieldData","Field Dataset Merged","CSV","FieldDataMerged Valid UTM.csv"))
list_Plot <- csv_Plot$Plot
# location of the directory where shapefiles used for validation are saved
vector_Dir_Sentinel <- list()
plot_res_Sentinel <- c(50,100,150,300)
for (x in 1:length(plot_res_Sentinel)){
  temp_Vector_Dir_Sentinel <- file.path(path_abs,"FieldData","Field Dataset Merged","Shapefiles Buffered",paste0("FieldDataMerged Valid Buffer ",plot_res_Sentinel[[x]],"m UTM.shp"))
  vector_Dir_Sentinel <- append(vector_Dir_Sentinel,list(temp_Vector_Dir_Sentinel))
  rm(temp_Vector_Dir_Sentinel)
  }
print("----READ FIELD PLOT SHP DONE----")

### Calculate Plot Biodiversity ###
# Calculate those indices
biodiv_Indicators_Sentinel <- list()
for (x in 1:length(subdir_Sentinel)){
  biodiv_Indicators_temp_List <- list()
  for (y in 1:length(vector_Dir_Sentinel)){
    biodiv_Indicators_temp <- diversity_from_plots_nofunc(Raster_SpectralSpecies = list_spectral_Sentinel[[x]]$SpectralSpecies, 
                                                          Plots = vector_Dir_Sentinel[[y]],
                                                          nbclusters = nbClusters_Sentinel, 
                                                          Raster_Functional = list_PCASentinel[[x]]$PCA_Files, 
                                                          Selected_Features = pc_sel_Sentinel)
    biodiv_Indicators_temp_List <- append(biodiv_Indicators_temp_List, list(biodiv_Indicators_temp))
    rm(biodiv_Indicators_temp)
  }
  biodiv_Indicators_Sentinel <- append(biodiv_Indicators_Sentinel, list(biodiv_Indicators_temp_List))
  rm(biodiv_Indicators_temp_List)
}
# Save those indices to .csv files
print("----SAVING ALPHA BIODIVERSITY INDICES----")
for (x in 1:length(subdir_Sentinel)){
  for (y in 1:length(vector_Dir_Sentinel)){
    temp_Bio <- biodiv_Indicators_Sentinel[[x]][[y]]
    temp_Results <- data.frame(list_Plot, temp_Bio$Richness, temp_Bio$Fisher,
                          temp_Bio$Shannon, temp_Bio$Simpson)
    names(temp_Results)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson")
    write.table(temp_Results, file = file.path(outDir_Sentinel,name_Sentinel[[x]],paste0("AlphaDiversity",plot_res_Sentinel[[y]],"m.csv")),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
    rm(temp_Bio,temp_Results)
  }
}
print("----SAVED ALPHA BIODIVERSITY INDICES----")
print("----SAVING BETA BIODIVERSITY INDICES----")
for (x in 1:length(subdir_Sentinel)){
  for (y in 1:length(vector_Dir_Sentinel)){
    temp_BC <- biodiv_Indicators_Sentinel[[x]][[y]]$BCdiss
    write.table(temp_BC, file.path(outDir_Sentinel,name_Sentinel[[x]],paste0("BrayCurtis",plot_res_Sentinel[[y]],"m.csv")),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
    rm(temp_BC)
  }
}
print("----SAVED BETA BIODIVERSITY INDICES----")

rm(x,y)