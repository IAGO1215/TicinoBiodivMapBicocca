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
window_size <- list(5,10)
# computational parameters
nbCPU <- 4
MaxRAM <- 0.2
# number of clusters (spectral species)
nbClusters <- list(20,35,50)
# spectral filtering - for our case, we don't need any spectrals to be excluded
spectral_excluded <- NULL

### PCA ###

print("----PERFORMING PCA----")

list_PCA <- list()
for (x in 1:length(subdir_prodata)){
  PCA_Output <- biodivMapR::perform_PCA(Input_Image_File = path_raster[[x]],
                                        Output_Dir = file.path(outDir,"PCA"),
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

### Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

output_Path_Cluster <- list()
for (x in 1:length(nbClusters)){
  output_Path_Cluster <- append(output_Path_Cluster,file.path(outDir, paste0(nbClusters[[x]],'Clusters')))
}

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))
list_spectral <- list()
for (x in 1:length(subdir_prodata)){
  list_spectral_temp <- list()
  for (y in 1:length(nbClusters)){
    Kmeans_info <- map_spectral_species_py(Input_Image_File = path_raster[[x]],
                                                       Input_Mask_File = path_mask[[x]],
                                                       Output_Dir = output_Path_Cluster[[y]],
                                                       SpectralSpace_Output = list_PCA[[x]],
                                                       SelectedPCs = selected_PCs,
                                                       nbclusters = nbClusters[[y]],
                                                       nbCPU = nbCPU, 
                                                       MaxRAM = MaxRAM)
    print(paste0('Spectral species generated for ',path_raster[[x]],' cluster ',nbClusters[[y]]))
    list_spectral_temp <- append(list_spectral_temp, list(Kmeans_info))
  }
  list_spectral <- append(list_spectral, list(list_spectral_temp))
}

print("----SPECTRAL SPECIES COMPUTED----")

### Alpha Beta Diversity ###
print("----MAP ALPHA DIVERSITY----")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')

for (x in 1:length(subdir_prodata)){
  for (y in 1:length(nbClusters)){
    for (z in 1:length(window_size)){
      biodivMapR::map_alpha_div(Input_Image_File = path_raster[[x]],
                    Input_Mask_File = path_mask[[x]],
                    Output_Dir = output_Path_Cluster[[y]],
                    TypePCA = typePCA,
                    window_size = window_size[[z]],
                    nbCPU = nbCPU,
                    MaxRAM = MaxRAM,
                    Index_Alpha = Index_Alpha,
                    nbclusters = nbClusters[[y]])
      print(paste0('Alpha diversity map generated for ',path_raster[[x]],' Cluster ',output_Path_Cluster[[y]],' Window ',window_size[[z]]))
    }
  }
}

print("----MAP BETA DIVERSITY----")
for (x in 1:length(subdir_prodata)){
  for (y in 1:length(nbClusters)){
    for (z in 1:length(window_size)){
      biodivMapR::map_beta_div(Input_Image_File = path_raster[[x]],
                   Output_Dir = output_Path_Cluster[[y]],
                   TypePCA = typePCA,
                   window_size = window_size[[z]],
                   nbCPU = nbCPU,
                   MaxRAM = MaxRAM,
                   nbclusters = nbClusters[[y]])
      print(paste0('Beta diversity map generated for ',path_raster[[x]]))
    }
  }
}

print("----BETA DIVERSITY COMPUTED----")

### Field Splot Biodiversity ###

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