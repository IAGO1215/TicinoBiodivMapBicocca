### This R script is aimed at calculating the alpha and beta biodiversity indices from the PRISMA raster data ###
### The workflow is: import rasters and mask rasters -> calculate PCA and select PCs -> calculate spectral species -> calculate alpha and beta ###

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
dir_PRISMAMerged <- file.path(path_abs, 'PRISMA Raster Raw', 'Merged')
# Get the names of all the files
name_PRISMAMerged <- unique(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir_PRISMAMerged, full.name = FALSE, recursive = FALSE))))
# The mask folder
dir_PRISMAMask <- file.path(path_abs, 'PRISMA Raster Mask')

# Get the paths to all the raster files! 
path_PRISMAMerged <- list()
for (x in 1:length(name_PRISMAMerged)) {
  temp <- file.path(dir_PRISMAMerged, name_PRISMAMerged[[x]])
  print(temp)
  path_PRISMAMerged <- append(path_PRISMAMerged, temp)
  rm(temp)
}
# Mask paths (Not including the original masks which are irrelevant to the current project)
path_PRISMAMask <- list()
for (x in 1:length(name_PRISMAMerged)) {
  temp <- file.path(dir_PRISMAMask, paste0(name_PRISMAMerged[[x]]," Mask"))
  print(temp)
  path_PRISMAMask <- append(path_PRISMAMask, temp)
  rm(temp)
}
# Output directory (No need to create the subfolders on our own since biodiverR will create them automatically)
outDir_PRISMA <- file.path(path_abs,'ResultsPRISMA')

### Setting parameters ###

print("----SETTING PARAMETERS FOR BIODIVERSITY----")

# continuum removal
continuum_Removal <- TRUE
# Type of dimensionality reduction
typePCA <- 'SPCA'
# Automatically set to FALSE if TypePCA     = 'MNF'
filterPCA <- FALSE
# window size forcomputation of spectral diversity
window_size_PRISMA <- 5
# computational parameters
nbCPU <- 4
MaxRAM <- 0.2
# number of clusters (spectral species)
nbClusters_PRISMA <- 20
# spectral filtering - for our case, we don't need any spectrals to be excluded
spectral_excluded <- NULL

### PCA All 230 Bands ###

print("----PERFORMING PCA----")

list_PCAPRISMA <- list()
for (x in 1:length(name_PRISMAMerged)){
  PCA_Output_PRISMA <- biodivMapR::perform_PCA(Input_Image_File = path_PRISMAMerged[[x]],
                                        Output_Dir = outDir_PRISMA,
                                        Input_Mask_File = path_PRISMAMask[[x]],
                                        TypePCA = typePCA,
                                        NbPCs_To_Keep = 230,
                                        FilterPCA = filterPCA,
                                        nbCPU = nbCPU,
                                        MaxRAM = MaxRAM,
                                        Continuum_Removal = continuum_Removal, 
                                        Excluded_WL = spectral_excluded)
  print(paste0('PCA generated for ',path_PRISMAMerged[[x]]))
  list_PCAPRISMA <- append(list_PCAPRISMA, list(PCA_Output_PRISMA))
  rm(PCA_Output_PRISMA)
}

print("----PERFORMING PCA SUCCESSFULLY----")

#### Save eigenvalues! #### 

for (x in 1:length(name_PRISMAMerged)){
  temp_Summary <- list_PCAPRISMA[[x]]$PCA_model$sdev
  temp_Summary_df <- data.frame(temp_Summary)
  names(temp_Summary_df)  = c("Standard Deviation")
  write.table(temp_Summary, file = file.path(outDir_PRISMA,name_PRISMAMerged[[x]],"PCA R.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_Summary,temp_Summary_df)
}

### Select PCs Mannually ###

pc_sel_PRISMA_6 <- c(1,4,9)
pc_sel_PRISMA_9 <- c(1,2,4)
pc_sel_PRISMA_9_New <- c(1,3,5)
pc_sel_PRISMA <- list(pc_sel_PRISMA_6,pc_sel_PRISMA_9,pc_sel_PRISMA_9_New)

### Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))
list_spectral_PRISMA <- list()
for (x in 1:length(name_PRISMAMerged)){
  Kmeans_info_PRISMA <- map_spectral_species_py(Input_Image_File = path_PRISMAMerged[[x]],
                                         Input_Mask_File = path_PRISMAMask[[x]],
                                         Output_Dir = outDir_PRISMA,
                                         SpectralSpace_Output = list_PCAPRISMA[[x]],
                                         SelectedPCs = pc_sel_PRISMA[[x]],
                                         nbclusters = nbClusters_PRISMA,
                                         nbCPU = nbCPU, 
                                         MaxRAM = MaxRAM)
  print('Spectral species generated')
  list_spectral_PRISMA <- append(list_spectral_PRISMA, list(Kmeans_info_PRISMA))
  rm(Kmeans_info_PRISMA)
}

print("----SPECTRAL SPECIES COMPUTED----")

### Alpha Beta Diversity ###
print("----MAP ALPHA DIVERSITY----")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')

for (x in 1:length(name_PRISMAMerged)){
  biodivMapR::map_alpha_div(Input_Image_File = path_PRISMAMerged[[x]],
                            Input_Mask_File = path_PRISMAMask[[x]],
                            Output_Dir = outDir_PRISMA,
                            TypePCA = typePCA,
                            window_size = window_size_PRISMA,
                            nbCPU = nbCPU,
                            MaxRAM = MaxRAM,
                            Index_Alpha = Index_Alpha,
                            nbclusters = nbClusters_PRISMA)
  print('Alpha diversity map generated!')
}

print("----MAP BETA DIVERSITY----")
for (x in 1:length(name_PRISMAMerged)){
  biodivMapR::map_beta_div(Input_Image_File = path_PRISMAMerged[[x]],
                           Output_Dir = outDir_PRISMA,
                           TypePCA = typePCA,
                           window_size = window_size_PRISMA,
                           nbCPU = nbCPU,
                           MaxRAM = MaxRAM,
                           nbclusters = nbClusters_PRISMA)
  print('Beta diversity map generated!')
}

print("----BETA DIVERSITY COMPUTED----")

### Field Plot Biodiversity ###

print("----FIELD PLOT STEP----")
### 300m Buffer
# Get the name for all the plot
csv_Plot <- read.csv(file.path(path_abs,"FieldData","Field Dataset Merged","CSV","FieldDataMerged Valid UTM.csv"))
list_Plot <- csv_Plot$Plot
# location of the directory where shapefiles used for validation are saved
vector_Dir_PRISMA <- file.path(path_abs,"FieldData","Field Dataset Merged","Shapefiles Buffered","FieldDataMerged Valid Buffer 300m UTM.shp")
# list vector data (In our case, there is only one shapefile, so alternatively we can directly input its abs path)
# vector_Path <- biodivMapR::list_shp(vector_Dir)
# vector_Name <- tools::file_path_sans_ext(basename(vector_Path))
# Calculate those indices
Biodiv_Indicators_PRISMA <- list()
for (x in 1:length(name_PRISMAMerged)){
  Biodiv_Indicators_PRISMA_temp <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral_PRISMA[[x]]$SpectralSpecies, 
                                                               Plots = vector_Dir_PRISMA,
                                                               nbclusters = nbClusters_PRISMA, 
                                                               Raster_Functional = list_PCAPRISMA[[x]]$PCA_Files, 
                                                               Selected_Features = pc_sel_PRISMA[[x]])
  Biodiv_Indicators_PRISMA <- append(Biodiv_Indicators_PRISMA, list(Biodiv_Indicators_PRISMA_temp))
  rm(Biodiv_Indicators_PRISMA_temp)
}
# Save those indices to .csv files
for (x in 1:length(name_PRISMAMerged)){
  temp_Bio_PRISMA <- Biodiv_Indicators_PRISMA[[x]]
  temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                    temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson,
                                    temp_Bio_PRISMA$FunctionalDiversity$FRic,
                                    temp_Bio_PRISMA$FunctionalDiversity$FEve,
                                    temp_Bio_PRISMA$FunctionalDiversity$FDiv)
  names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
  write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,name_PRISMAMerged[[x]],"AlphaDiversity300m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_Bio_PRISMA,temp_Results_PRISMA)
}
for (x in 1:length(name_PRISMAMerged)){
  temp_BC_PRISMA <- Biodiv_Indicators_PRISMA[[x]]$BCdiss
  write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,name_PRISMAMerged[[x]],"BrayCurtis300m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_BC_PRISMA)
}
### 150m Buffer
vector_Dir_PRISMA5 <- file.path(path_abs,"FieldData","Field Dataset Merged","Shapefiles Buffered","FieldDataMerged Valid Buffer 150m UTM.shp")
# Calculate those indices
Biodiv_Indicators_PRISMA5 <- list()
for (x in 1:length(name_PRISMAMerged)){
  Biodiv_Indicators_PRISMA_temp <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral_PRISMA[[x]]$SpectralSpecies, 
                                                                    Plots = vector_Dir_PRISMA5,
                                                                    nbclusters = nbClusters_PRISMA, 
                                                                    Raster_Functional = list_PCAPRISMA[[x]]$PCA_Files, 
                                                                    Selected_Features = pc_sel_PRISMA[[x]])
  Biodiv_Indicators_PRISMA5 <- append(Biodiv_Indicators_PRISMA5, list(Biodiv_Indicators_PRISMA_temp))
  rm(Biodiv_Indicators_PRISMA_temp)
}
# Save those indices to .csv files
for (x in 1:length(name_PRISMAMerged)){
  temp_Bio_PRISMA <- Biodiv_Indicators_PRISMA5[[x]]
  temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                    temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson,
                                    temp_Bio_PRISMA$FunctionalDiversity$FRic,
                                    temp_Bio_PRISMA$FunctionalDiversity$FEve,
                                    temp_Bio_PRISMA$FunctionalDiversity$FDiv)
  names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
  write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,name_PRISMAMerged[[x]],"AlphaDiversity150m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_Bio_PRISMA,temp_Results_PRISMA)
}
for (x in 1:length(name_PRISMAMerged)){
  temp_BC_PRISMA <- Biodiv_Indicators_PRISMA5[[x]]$BCdiss
  write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,name_PRISMAMerged[[x]],"BrayCurtis150m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_BC_PRISMA)
}

### Field Splot Biodiversity with 15 PCs Spectral Species###

print("----FIELD PLOT STEP----")
### 300m Buffer

# Calculate those indices
Biodiv_Indicators_PRISMA10_15 <- list()
for (x in 1:length(name_PRISMAMerged)){
  Biodiv_Indicators_PRISMA_temp <- diversity_from_plots_nofunc(Raster_SpectralSpecies = list_spectral_PRISMA_15[[x]]$SpectralSpecies, 
                                                                    Plots = vector_Dir_PRISMA,
                                                                    nbclusters = nbClusters_PRISMA, 
                                                                    Raster_Functional = list_PCAPRISMA[[x]]$PCA_Files, 
                                                                    Selected_Features = pc_sel_PRISMA15)
  Biodiv_Indicators_PRISMA10_15 <- append(Biodiv_Indicators_PRISMA10_15, list(Biodiv_Indicators_PRISMA_temp))
}
# Save those indices to .csv files
for (x in 1:length(name_PRISMAMerged)){
  temp_Bio_PRISMA <- Biodiv_Indicators_PRISMA10_15[[x]]
  temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                    temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson)
  names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson")
  write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,"15",name_PRISMAMerged[[x]],"AlphaDiversity10_15PCs.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
}
for (x in 1:length(name_PRISMAMerged)){
  temp_BC_PRISMA <- Biodiv_Indicators_PRISMA10_15[[x]]$BCdiss
  write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,"15",name_PRISMAMerged[[x]],"BrayCurtis10_15PCs.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
}
# Calculate those indices
Biodiv_Indicators_PRISMA5_15 <- list()
for (x in 1:length(name_PRISMAMerged)){
  Biodiv_Indicators_PRISMA_temp <- diversity_from_plots_nofunc(Raster_SpectralSpecies = list_spectral_PRISMA_15[[x]]$SpectralSpecies, 
                                                               Plots = vector_Dir_PRISMA5,
                                                               nbclusters = nbClusters_PRISMA, 
                                                               Raster_Functional = list_PCAPRISMA[[x]]$PCA_Files, 
                                                               Selected_Features = pc_sel_PRISMA15)
  Biodiv_Indicators_PRISMA5_15 <- append(Biodiv_Indicators_PRISMA5_15, list(Biodiv_Indicators_PRISMA_temp))
}
# Save those indices to .csv files
for (x in 1:length(name_PRISMAMerged)){
  temp_Bio_PRISMA <- Biodiv_Indicators_PRISMA5_15[[x]]
  temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                    temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson)
  names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson")
  write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,"15",name_PRISMAMerged[[x]],"AlphaDiversity5_15PCs.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
}
for (x in 1:length(name_PRISMAMerged)){
  temp_BC_PRISMA <- Biodiv_Indicators_PRISMA5_15[[x]]$BCdiss
  write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,"15",name_PRISMAMerged[[x]],"BrayCurtis5_15PCs.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
}