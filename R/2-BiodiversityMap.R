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

### Field Plot Biodiversity ###

print("----FIELD PLOT STEP----")
# Get the name for all the plot
csv_Plot <- read.csv(file.path(path_abs,"FieldData","Field Dataset Merged","CSV","FieldDataMerged Valid UTM.csv"))
list_Plot <- csv_Plot$Plot
# location of the directory where shapefiles used for validation are saved
vector_Dir <- file.path(path_abs,"FieldData","Field Dataset Merged","Shapefiles Buffered","FieldDataMerged Valid Buffer 50m UTM.shp")
# list vector data (In our case, there is only one shapefile, so alternatively we can directly input its abs path)
# vector_Path <- biodivMapR::list_shp(vector_Dir)
# vector_Name <- tools::file_path_sans_ext(basename(vector_Path))
# Calculate those indices
Biodiv_Indicators <- list()
for (x in 1:length(list_spectral)){
  Biodiv_Indicators_temp_List <- list()
  for (y in 1:length(nbClusters)){
    Biodiv_Indicators_temp <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral[[x]][[y]]$SpectralSpecies, 
                                                          Plots = vector_Dir,
                                                          nbclusters = nbClusters[[y]], 
                                                          Raster_Functional = list_PCA[[x]]$PCA_Files, 
                                                          Selected_Features = selected_PCs)
    Biodiv_Indicators_temp_List <- append(Biodiv_Indicators_temp_List, list(Biodiv_Indicators_temp))
  }
  Biodiv_Indicators <- append(Biodiv_Indicators, list(Biodiv_Indicators_temp_List))
}
# Save those indices to .csv files
for (x in 1:length(list_spectral)){
  for (y in 1:length(nbClusters)){
    temp_Bio <- Biodiv_Indicators[[x]][[y]]
    temp_Results <- data.frame(list_Plot, temp_Bio$Richness, temp_Bio$Fisher,
                          temp_Bio$Shannon, temp_Bio$Simpson,
                          temp_Bio$FunctionalDiversity$FRic,
                          temp_Bio$FunctionalDiversity$FEve,
                          temp_Bio$FunctionalDiversity$FDiv)
    names(temp_Results)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
    write.table(temp_Results, file = file.path(output_Path_Cluster[[y]],basename(path_raster[[x]]),"AlphaDiversity.csv"),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  }
}
for (x in 1:length(list_spectral)){
  for (y in 1:length(nbClusters)){
    temp_BC <- Biodiv_Indicators[[x]][[y]]$BCdiss
    write.table(temp_BC, file.path(output_Path_Cluster[[y]],basename(path_raster[[x]]),"BrayCurtis.csv"),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  }
}
# 10x10 plot
vector_Dir_10 <- file.path(path_abs,"FieldData","Field Dataset Merged","Shapefiles Buffered","FieldDataMerged Valid Buffer 100m UTM.shp")
Biodiv_Indicators_10 <- list()
for (x in 1:length(list_spectral)){
  Biodiv_Indicators_temp_List <- list()
  for (y in 1:length(nbClusters)){
    Biodiv_Indicators_temp <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral[[x]][[y]]$SpectralSpecies, 
                                                               Plots = vector_Dir_10,
                                                               nbclusters = nbClusters[[y]], 
                                                               Raster_Functional = list_PCA[[x]]$PCA_Files, 
                                                               Selected_Features = selected_PCs)
    Biodiv_Indicators_temp_List <- append(Biodiv_Indicators_temp_List, list(Biodiv_Indicators_temp))
  }
  Biodiv_Indicators_10 <- append(Biodiv_Indicators_10, list(Biodiv_Indicators_temp_List))
}
# Save those indices to .csv files
for (x in 1:length(list_spectral)){
  for (y in 1:length(nbClusters)){
    temp_Bio <- Biodiv_Indicators_10[[x]][[y]]
    temp_Results <- data.frame(list_Plot, temp_Bio$Richness, temp_Bio$Fisher,
                               temp_Bio$Shannon, temp_Bio$Simpson,
                               temp_Bio$FunctionalDiversity$FRic,
                               temp_Bio$FunctionalDiversity$FEve,
                               temp_Bio$FunctionalDiversity$FDiv)
    names(temp_Results)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
    write.table(temp_Results, file = file.path(output_Path_Cluster[[y]],basename(path_raster[[x]]),"AlphaDiversity10.csv"),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  }
}
for (x in 1:length(list_spectral)){
  for (y in 1:length(nbClusters)){
    temp_BC <- Biodiv_Indicators_10[[x]][[y]]$BCdiss
    write.table(temp_BC, file.path(output_Path_Cluster[[y]],basename(path_raster[[x]]),"BrayCurtis10.csv"),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  }
}