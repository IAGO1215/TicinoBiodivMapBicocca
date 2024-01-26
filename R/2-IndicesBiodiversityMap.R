# Import the R package
library(biodivMapR)
library(progress)
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
path_Sentinel_RasterIndices <- list()
for (x in 1:length(subdir_Sentinel)) {
  temp <- file.path(subdir_Sentinel[x], paste0(name_Sentinel[x],'Cropped_StackedIndices_20VI'))
  print(temp)
  path_Sentinel_RasterIndices <- append(path_Sentinel_RasterIndices, temp)
  rm(temp)
}
### PCA ###

print("----PERFORMING PCA----")

list_PCASentinelIndices <- list()
for (x in 1:length(subdir_Sentinel)){
  temp_PCA <- biodivMapR::perform_PCA(Input_Image_File = path_Sentinel_RasterIndices[[x]],
                                      Output_Dir = outDir_Sentinel,
                                      Input_Mask_File = path_Sentinel_Mask[[x]],
                                      TypePCA = typePCA,
                                      FilterPCA = filterPCA,
                                      nbCPU = nbCPU,
                                      MaxRAM = MaxRAM,
                                      Continuum_Removal = FALSE, 
                                      Excluded_WL = spectral_excluded)
  print(paste0('PCA generated for ',path_Sentinel_RasterIndices[[x]]))
  list_PCASentinelIndices <- append(list_PCASentinelIndices, list(temp_PCA))
  rm(temp_PCA)
}

print("----PERFORMING PCA SUCCESSFULLY----")

print("----MANUALLY SELECTING PCs----")
pc_sel_SentinelIndices <- c(1,2,3,4,5,6,7,8,9,10)
### Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))

list_spectral_SentinelIndices <- list()
for (x in 1:length(subdir_Sentinel)){
  Kmeans_info <- map_spectral_species_py(Input_Image_File = path_Sentinel_RasterIndices[[x]],
                                         Input_Mask_File = path_Sentinel_Mask[[x]],
                                         Output_Dir = outDir_Sentinel,
                                         SpectralSpace_Output = list_PCASentinelIndices[[x]],
                                         SelectedPCs = pc_sel_SentinelIndices,
                                         nbclusters = nbClusters_Sentinel,
                                         nbCPU = nbCPU, 
                                         MaxRAM = MaxRAM)
  print(paste0('Spectral species generated for ',path_Sentinel_RasterIndices[[x]]))
  list_spectral_SentinelIndices <- append(list_spectral_SentinelIndices, list(Kmeans_info))
  rm(Kmeans_info)
}

print("----SPECTRAL SPECIES COMPUTED----")

### Alpha Beta Diversity ###
print("----MAP ALPHA DIVERSITY----")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')

for (x in 1:length(subdir_Sentinel)){
  biodivMapR::map_alpha_div(Input_Image_File = path_Sentinel_RasterIndices[[x]],
                            Input_Mask_File = path_Sentinel_Mask[[x]],
                            Output_Dir = outDir_Sentinel,
                            TypePCA = typePCA,
                            window_size = window_size_Sentinel,
                            nbCPU = nbCPU,
                            MaxRAM = MaxRAM,
                            Index_Alpha = Index_Alpha,
                            nbclusters = nbClusters_Sentinel)
  print(paste0('Alpha diversity map generated for ',path_Sentinel_RasterIndices[[x]]))
}

print("----MAP BETA DIVERSITY----")
for (x in 1:length(subdir_Sentinel)){
  biodivMapR::map_beta_div(Input_Image_File = path_Sentinel_RasterIndices[[x]],
                           Output_Dir = outDir_Sentinel,
                           TypePCA = typePCA,
                           window_size = window_size_Sentinel,
                           nbCPU = nbCPU,
                           MaxRAM = MaxRAM,
                           nbclusters = nbClusters_Sentinel)
  print(paste0('Beta diversity map generated for ',path_Sentinel_RasterIndices[[x]]))
}

print("----BETA DIVERSITY COMPUTED----")

### Calculate Plot Biodiversity ###
# Calculate those indices
biodiv_Indicators_SentinelIndices <- list()
for (x in 1:length(subdir_Sentinel)){
  biodiv_Indicators_temp_List <- list()
  for (y in 1:length(vector_Dir_Sentinel)){
    biodiv_Indicators_temp <- diversity_from_plots_nofunc(Raster_SpectralSpecies = list_spectral_SentinelIndices[[x]]$SpectralSpecies, 
                                                          Plots = vector_Dir_Sentinel[[y]],
                                                          nbclusters = nbClusters_Sentinel, 
                                                          Raster_Functional = list_PCASentinelIndices[[x]]$PCA_Files, 
                                                          Selected_Features = pc_sel_SentinelIndices)
    biodiv_Indicators_temp_List <- append(biodiv_Indicators_temp_List, list(biodiv_Indicators_temp))
    rm(biodiv_Indicators_temp)
  }
  biodiv_Indicators_SentinelIndices <- append(biodiv_Indicators_SentinelIndices, list(biodiv_Indicators_temp_List))
  rm(biodiv_Indicators_temp_List)
}
# Save those indices to .csv files
print("----SAVING ALPHA BIODIVERSITY INDICES----")
for (x in 1:length(subdir_Sentinel)){
  for (y in 1:length(vector_Dir_Sentinel)){
    temp_Bio <- biodiv_Indicators_SentinelIndices[[x]][[y]]
    temp_Results <- data.frame(list_Plot, temp_Bio$Richness, temp_Bio$Fisher,
                               temp_Bio$Shannon, temp_Bio$Simpson)
    names(temp_Results)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson")
    write.table(temp_Results, file = file.path(outDir_Sentinel,paste0(name_Sentinel[x],'Cropped_StackedIndices_20VI'),paste0("AlphaDiversity",plot_res_Sentinel[[y]],"m.csv")),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
    rm(temp_Bio,temp_Results)
  }
}
print("----SAVED ALPHA BIODIVERSITY INDICES----")
print("----SAVING BETA BIODIVERSITY INDICES----")
for (x in 1:length(subdir_Sentinel)){
  for (y in 1:length(vector_Dir_Sentinel)){
    temp_BC <- biodiv_Indicators_SentinelIndices[[x]][[y]]$BCdiss
    write.table(temp_BC, file.path(outDir_Sentinel,paste0(name_Sentinel[x],'Cropped_StackedIndices_20VI'),paste0("BrayCurtis",plot_res_Sentinel[[y]],"m.csv")),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
    rm(temp_BC)
  }
}
print("----SAVED BETA BIODIVERSITY INDICES----")

rm(x,y)
