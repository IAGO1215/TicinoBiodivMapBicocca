# Import the R package
library(biodivMapR)
library(spinR)
library(stars)
# Set absolute working directory
path_abs <- "C:/Users/m1865/Desktop/Ticino"
setwd(path_abs)

# Calculate indices and stack them! 
path_PRISMAIndices <- list()
for (x in 1:length(name_PRISMAMerged)){
  ### Calculate spectral indices ### 
  ImBrick <- raster::brick(path_PRISMAMerged[[x]])
  HDR <- read_ENVI_header(get_HDR_name(path_PRISMAMerged[[x]]))
  SensorBands <- HDR$wavelength
  # compute the spectral indices from raster data, using spectral bands as 
  # close as possible from Sentinel-2 spectral bands used for the spectral indices
  all_Indices <- list("ARI1","ARI2","ARVI","BAI","BAIS2","CCCI","CHL_RE","CRI1","CRI2","EVI","EVI2","GRVI1","GNDVI","IRECI",
                      "LAI_SAVI","MCARI","mNDVI705","MSAVI2","MSI","mSR705","MTCI","nBR_RAW", 
                       "NDI_45","NDII","NDSI","NDVI","NDVI_G","NDVI705","NDWI1","NDWI2","PSRI",
                      "PSRI_NIR","RE_NDVI","RE_NDWI","S2REP","SAVI","SIPI","SR","CR_SWIR","CR_RE" )
  Spectral_Indices <- spinR::compute_S2SI_Raster(Refl = ImBrick, 
                                                 SensorBands = SensorBands, 
                                                 Sel_Indices = all_Indices, 
                                                 StackOut = T, 
                                                 ReflFactor = 10000)
  # Mask
  # initialize mask
  Mask <- 0*Spectral_Indices$SpectralIndices$mNDVI705+1
  # remove outliers from spectral indices
  for (idx in all_Indices){
    rast <- Spectral_Indices$SpectralIndices[[idx]]
    IQRminmax <- biodivMapR::IQR_outliers(DistVal = raster::values(rast),weightIRQ = 3)
    Mask[raster::values(rast)<IQRminmax[1] | raster::values(rast)>IQRminmax[2]] <- NA
  }
  # apply mask on stack and convert as stars object
  StarsObj <- st_as_stars(Spectral_Indices$SpectralIndices*Mask)
  # Save mask
  # define path where to store spectral indices (same as root path for biodivMapR output directory)
  NameStack <- paste(name_PRISMAMerged[[x]],'_StackedIndices',sep = '')
  PathIndices <- file.path(outDir_PRISMA,NameStack)
  dir.create(PathIndices,recursive = T,showWarnings = F)
  # save mask in raster file
  Input_Mask_File_temp <- file.path(PathIndices,'Mask')
  stars::write_stars(st_as_stars(Mask), dsn=Input_Mask_File_temp, driver =  "ENVI",type='Byte')
  # save Stack of spectral indices in raster file with spectral indices defining band names
  Path_StackIndices <- file.path(PathIndices,NameStack)
  path_PRISMAIndices <- append (path_PRISMAIndices, Path_StackIndices)
  biodivMapR::write_StarsStack(StarsObj = StarsObj, 
                               dsn = Path_StackIndices, 
                               BandNames = all_Indices, 
                               datatype='Float32')
  rm(ImBrick,HDR,SensorBands,Spectral_Indices,Mask,rast,IQRminmax,StarsObj,NameStack,PathIndices,Input_Mask_File_temp,Path_StackIndices)
  }

### Use PCA Alternative as Input for Spectral Species ###
print("----COMPUTING SPECTRAL SPECIES----")

library(raster)
source(file.path(path_abs,'R','Updated_Functions.R'))
list_spectral_PRISMAIndices <- list()

for (x in 1:length(path_PRISMAIndices)){
  SpectralSpace_Output_Indices <- list('PCA_Files' = path_PRISMAIndices[[x]], 
                                       'TypePCA' = 'noPCA')
  Kmeans_info_PRISMA <- map_spectral_species_py(Input_Image_File = path_PRISMAIndices[[x]],
                                                Input_Mask_File = path_PRISMAMask[[x]],
                                                Output_Dir = outDir_PRISMA,
                                                SpectralSpace_Output = SpectralSpace_Output_Indices,
                                                SelectedPCs = c(6,16,26),
                                                nbclusters = 20,
                                                nbCPU = nbCPU, 
                                                MaxRAM = MaxRAM)
  print('Spectral species generated')
  list_spectral_PRISMAIndices <- append(list_spectral_PRISMAIndices, list(Kmeans_info_PRISMA))
  print(paste0(path_PRISMAIndices[[x]]," Saved!"))
  rm(SpectralSpace_Output_Indices,Kmeans_info_PRISMA)
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
                            TypePCA = 'noPCA',
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
                           TypePCA = 'noPCA',
                           window_size = window_size_PRISMA,
                           nbCPU = nbCPU,
                           MaxRAM = MaxRAM,
                           nbclusters = nbClusters_PRISMA)
  print('Beta diversity map generated!')
}

print("----BETA DIVERSITY COMPUTED----")

### Field Splot Biodiversity ###

print("----FIELD PLOT STEP----")
### 300m Buffer ###

# Calculate those indices
Biodiv_Indicators_PRISMAIN_300 <- list()
for (x in 1:length(path_PRISMAIndices)){
  Biodiv_Indicators_PRISMA_temp <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral_PRISMAIndices[[x]]$SpectralSpecies, 
                                                                    Plots = vector_Dir_PRISMA,
                                                                    nbclusters = nbClusters_PRISMA, 
                                                                    Raster_Functional = path_PRISMAIndices[[x]], 
                                                                    Selected_Features = c(6,16,26))
  Biodiv_Indicators_PRISMAIN_300 <- append(Biodiv_Indicators_PRISMAIN_300, list(Biodiv_Indicators_PRISMA_temp))
  rm(Biodiv_Indicators_PRISMA_temp)
}
# Save those indices to .csv files
for (x in 1:length(path_PRISMAIndices)){
  temp_Bio_PRISMA <- Biodiv_Indicators_PRISMAIN_300[[x]]
  temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                    temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson,
                                    temp_Bio_PRISMA$FunctionalDiversity$FRic,
                                    temp_Bio_PRISMA$FunctionalDiversity$FEve,
                                    temp_Bio_PRISMA$FunctionalDiversity$FDiv)
  names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
  write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,paste(name_PRISMAMerged[[x]],'_StackedIndices',sep = ''),"AlphaDiversity300m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_Bio_PRISMA,temp_Results_PRISMA)
}
for (x in 1:length(path_PRISMAIndices)){
  temp_BC_PRISMA <- Biodiv_Indicators_PRISMAIN_300[[x]]$BCdiss
  write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,paste(name_PRISMAMerged[[x]],'_StackedIndices',sep = ''),"BrayCurtis300m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_BC_PRISMA)
}

### 150m Buffer ###
# Calculate those indices
Biodiv_Indicators_PRISMAIN_150 <- list()
for (x in 1:length(path_PRISMAIndices)){
  Biodiv_Indicators_PRISMA_temp <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = list_spectral_PRISMAIndices[[x]]$SpectralSpecies, 
                                                                    Plots = vector_Dir_PRISMA5,
                                                                    nbclusters = nbClusters_PRISMA, 
                                                                    Raster_Functional = path_PRISMAIndices[[x]], 
                                                                    Selected_Features = c(6,16,26))
  Biodiv_Indicators_PRISMAIN_150 <- append(Biodiv_Indicators_PRISMAIN_150, list(Biodiv_Indicators_PRISMA_temp))
  rm(Biodiv_Indicators_PRISMA_temp)
}
# Save those indices to .csv files
for (x in 1:length(path_PRISMAIndices)){
  temp_Bio_PRISMA <- Biodiv_Indicators_PRISMAIN_150[[x]]
  temp_Results_PRISMA <- data.frame(list_Plot, temp_Bio_PRISMA$Richness, temp_Bio_PRISMA$Fisher,
                                    temp_Bio_PRISMA$Shannon, temp_Bio_PRISMA$Simpson,
                                    temp_Bio_PRISMA$FunctionalDiversity$FRic,
                                    temp_Bio_PRISMA$FunctionalDiversity$FEve,
                                    temp_Bio_PRISMA$FunctionalDiversity$FDiv)
  names(temp_Results_PRISMA)  = c("Plot","Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
  write.table(temp_Results_PRISMA, file = file.path(outDir_PRISMA,paste(name_PRISMAMerged[[x]],'_StackedIndices',sep = ''),"AlphaDiversity150m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_Bio_PRISMA,temp_Results_PRISMA)
}
for (x in 1:length(path_PRISMAIndices)){
  temp_BC_PRISMA <- Biodiv_Indicators_PRISMAIN_150[[x]]$BCdiss
  write.table(temp_BC_PRISMA, file.path(outDir_PRISMA,paste(name_PRISMAMerged[[x]],'_StackedIndices',sep = ''),"BrayCurtis150m.csv"),
              sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
  rm(temp_BC_PRISMA)
}