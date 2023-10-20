### This R script is aimed at converting original Sentinel-2 data (.SAFE) to bands-stacked raster file (ENVI header) ###
### All the .SAFE should be put into separate subfolders (can be renamed) under "RawData" ###
### All the output ENVI header raster file will be saved into corresponding automatically-created subfolders under "ProcessedData" ### 

# Clean environment
rm(list=ls(all=TRUE));gc()
# Import the R package
library(preprocS2)
# Set absolute working directory
path_abs <- "C:/Users/m1865/Desktop/Ticino"
setwd(path_abs)

### INPUT AND OUTPUT FOLDERS ###

print("----PREPARING FOLDERS FOR INPUT AND OUTPUT----")

# The raw data folder
dir_rawdata <- file.path(path_abs, 'RawData')
# Get all the subfolders in our "rawdata" folder
subdir_rawdata <- list.dirs(dir_rawdata, recursive = FALSE)
# Get the names of all the subfolders
name_rawdata <- list.dirs(dir_rawdata, full.name = FALSE, recursive = FALSE)
# The output data folder
dir_result <- file.path(path_abs, 'ProcessedData')
# Create the subfolders where the output reflactance file will be saved separately
subdir_result <- list()
for (x in name_rawdata) {
  subdir_result <- append(subdir_result, file.path(dir_result,x))
}

### PREPROCESSING

#### STACK INDIVIDUAL BANDS

print("----STACKING INDIVIDUAL BANDS----")

# Define resolution (Sentinel-2 resolution 10m x 10m)
resolution <- 10
# Define source of data
S2source <- 'SAFE'

# Stack individual bands for all the rawdata in different folders
list_stacked <- list()

for (x in subdir_rawdata) {
   print(paste0("----Stacking individual bands initiated for ", basename(x), "----"))  
   temp_stack <- preprocS2::extract_from_S2_L2A(Path_dir_S2 = x,
                                                S2source = S2source,
                                                resolution = resolution)
   print(paste0("----Stacking individual bands completed for ", basename(x), "----"))
   list_stacked <- append(list_stacked, list(temp_stack))
}

names(list_stacked) <- name_rawdata

print("----STACKING BANDS COMPLETED----")

#### CREATE REFLACTANCE DATA

print("----CREATING REFLACTANCE RASTER FILE----")

counter <- 0

for (x in list_stacked) {
  counter <- counter + 1
  print(paste0("Counter now equals to ", counter))
  # Prepare function input parameters
  s2_stars <- x$S2_Stack
  refl_path <- file.path(subdir_result[[counter]],basename(subdir_result[[counter]]))
  dir.create(subdir_result[[counter]])
  tile_S2 <- get_tile(x$S2_Bands$GRANULE)
  dateAcq_S2 <- get_date(x$S2_Bands$GRANULE)
  mtd <- x$S2_Bands$metadata
  mtd_msi <- x$S2_Bands$metadata_MSI
  
  # Call the function to produce reflactance data
  print(paste0("----Producing reflactance file initiated for ", basename(subdir_result[[counter]]), "----"))
  preprocS2::save_reflectance_s2(S2_stars = s2_stars,
                                 Refl_path = refl_path,
                                 S2Sat = NULL,
                                 tile_S2 = tile_S2,
                                 dateAcq_S2 = dateAcq_S2,
                                 Format = 'ENVI',
                                 datatype = 'Int16',
                                 MTD = mtd,
                                 MTD_MSI = mtd_msi,
                                 s2mission = '2A')
  print(paste0("----Producing reflactance file completed for ", basename(subdir_result[[counter]]), "----"))
}

print("----REFLACTANCE RASTER FILE CREATED SUCCESSFULLY----")

print(paste0("----The output reflactance file is saved at ",dir_result))
