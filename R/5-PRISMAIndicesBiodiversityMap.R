#### This R script calculates the PCA, alpha and beta biodiversity indices from the stacked 20 vegetation indices raster files.  ####

# Import the R package
library(biodivMapR)
# Set absolute working directory
path_abs <- "C:/Users/m1865/Desktop/Ticino"
setwd(path_abs)

### Reading files ###

print("----LOADING INPUT SATELLITE IMAGE FILES----")

# Importing raster files
# The raw data folder
dir_PRISMAIndices <- file.path(path_abs, 'PRISMA Classification Indices', 'Sel')
# Get the names of all the files
name_PRISMAIndices <- unique(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir_PRISMAIndices, full.name = FALSE, recursive = FALSE))))
name_PRISMAIndices_Whole <- character()
for (x in 1:length(name_PRISMAIndices)){
  if (!grepl("Cropped",name_PRISMAIndices[[x]])){
    name_PRISMAIndices_Whole <- c(name_PRISMAIndices_Whole,name_PRISMAIndices[[x]])
  }
}

# Get the paths to all the raster files! 
path_PRISMAIndices <- list()
for (x in 1:length(name_PRISMAIndices_Whole)) {
  temp <- file.path(dir_PRISMAIndices, name_PRISMAIndices_Whole[[x]])
  print(temp)
  path_PRISMAIndices <- append(path_PRISMAIndices, temp)
  rm(temp)
}

### PCA All Bands ###

print("----PERFORMING PCA----")

list_PCAPRISMAIndices <- list()
for (x in 1:length(path_PRISMAIndices)){
  PCA_Output_PRISMA <- biodivMapR::perform_PCA(Input_Image_File = path_PRISMAIndices[[x]],
                                               Output_Dir = outDir_PRISMA,
                                               Input_Mask_File = path_PRISMAMask[[x]],
                                               TypePCA = typePCA,
                                               FilterPCA = filterPCA,
                                               nbCPU = nbCPU,
                                               MaxRAM = MaxRAM,
                                               Continuum_Removal = FALSE, 
                                               Excluded_WL = spectral_excluded)
  print(paste0('PCA generated for ',name_PRISMAIndices_Whole[[x]]))
  list_PCAPRISMAIndices <- append(list_PCAPRISMAIndices, list(PCA_Output_PRISMA))
  rm(PCA_Output_PRISMA)
}

print("----PERFORMING PCA SUCCESSFULLY----")

#### Save eigenvalues! #### 

for (x in 1:length(path_PRISMAIndices)){
  temp_Summary <- list_PCAPRISMAIndices[[x]]$PCA_model$sdev
  temp_Summary_df <- data.frame(temp_Summary)
  names(temp_Summary_df)  = c("Standard Deviation")
  write.table(temp_Summary, file = file.path(outDir_PRISMA,name_PRISMAIndices_Whole[[x]],"PCA R.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_Summary,temp_Summary_df)
}

#### Select Bands ####
pc_sel_PRISMAIndices <- c(1,2,3,4,5,6,7,8,9,10)

### Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))
list_spectral_PRISMAIndices <- list()
for (x in 1:length(path_PRISMAIndices)){
  Kmeans_info_PRISMA <- map_spectral_species_py(Input_Image_File = path_PRISMAIndices[[x]],
                                                Input_Mask_File = path_PRISMAMask[[x]],
                                                Output_Dir = outDir_PRISMA,
                                                SpectralSpace_Output = list_PCAPRISMAIndices[[x]],
                                                SelectedPCs = pc_sel_PRISMAIndices,
                                                nbclusters = nbClusters_PRISMA,
                                                nbCPU = nbCPU, 
                                                MaxRAM = MaxRAM)
  print('Spectral species generated')
  list_spectral_PRISMAIndices <- append(list_spectral_PRISMAIndices, list(Kmeans_info_PRISMA))
  rm(Kmeans_info_PRISMA)
}

print("----SPECTRAL SPECIES COMPUTED----")

### Alpha Beta Diversity ###
print("----MAP ALPHA DIVERSITY----")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')

for (x in 1:length(path_PRISMAIndices)){
  biodivMapR::map_alpha_div(Input_Image_File = path_PRISMAIndices[[x]],
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
for (x in 1:length(path_PRISMAIndices)){
  biodivMapR::map_beta_div(Input_Image_File = path_PRISMAIndices[[x]],
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
# Get the name for all the plot
csv_Plot <- read.csv(file.path(path_abs,"FieldData","Field Dataset Merged","CSV","FieldDataMerged Valid UTM.csv"))
list_Plot <- csv_Plot$Plot
# location of the directory where shapefiles used for validation are saved
vector_Dir_PRISMA <- list()
plot_res_PRISMA <- c(150,300)
for (x in 1:length(plot_res_PRISMA)){
  temp_Vector_Dir_PRISMA <- file.path(path_abs,"FieldData","Field Dataset Merged","Shapefiles Buffered",paste0("FieldDataMerged Valid Buffer ",plot_res_PRISMA[[x]],"m UTM.shp"))
  vector_Dir_PRISMA <- append(vector_Dir_PRISMA,list(temp_Vector_Dir_PRISMA))
  rm(temp_Vector_Dir_PRISMA)
}
# list vector data (In our case, there is only one shapefile, so alternatively we can directly input its abs path)
# vector_Path <- biodivMapR::list_shp(vector_Dir)
# vector_Name <- tools::file_path_sans_ext(basename(vector_Path))
# Calculate those indices
biodiv_Indicators_PRISMAIndices <- list()
for (x in 1:length(path_PRISMAIndices)){
  biodiv_Indicators_temp_List <- list()
  for (y in 1:length(vector_Dir_PRISMA)){
    biodiv_Indicators_temp <- diversity_from_plots_nofunc(Raster_SpectralSpecies = list_spectral_PRISMAIndices[[x]]$SpectralSpecies, 
                                                          Plots = vector_Dir_PRISMA[[y]],
                                                          nbclusters = nbClusters_PRISMA, 
                                                          Raster_Functional = list_PCAPRISMAIndices[[x]]$PCA_Files, 
                                                          Selected_Features = pc_sel_PRISMAIndices)
    biodiv_Indicators_temp_List <- append(biodiv_Indicators_temp_List, list(biodiv_Indicators_temp))
    rm(biodiv_Indicators_temp)
  }
  biodiv_Indicators_PRISMAIndices <- append(biodiv_Indicators_PRISMAIndices, list(biodiv_Indicators_temp_List))
  rm(biodiv_Indicators_temp_List)
}
# Save those indices to .csv files
for (x in 1:length(path_PRISMAIndices)){
  for (y in 1:length(vector_Dir_PRISMA)){
    temp_Bio_PRISMA <- biodiv_Indicators_PRISMAIndices[[x]][[y]]
    temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                      temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson)
    names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson")
    write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,name_PRISMAIndices_Whole[[x]],paste0("AlphaDiversity", plot_res_PRISMA[[y]], "m.csv")),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
    rm(temp_Bio_PRISMA,temp_Results_PRISMA)
  }
}
for (x in 1:length(path_PRISMAIndices)){
  for (y in 1:length(vector_Dir_PRISMA)){
    temp_BC_PRISMA <- biodiv_Indicators_PRISMAIndices[[x]][[y]]$BCdiss
    write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,name_PRISMAIndices_Whole[[x]],paste0("BrayCurtis",plot_res_PRISMA[[y]],"m.csv")),
                sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
    rm(temp_BC_PRISMA)
  }
}