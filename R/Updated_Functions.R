library(raster)
library(biodivMapR)
library(tools)
library(progress)

map_spectral_species_py <- function (Input_Image_File, Output_Dir, SpectralSpace_Output, 
          Input_Mask_File = FALSE, nbclusters = 50, nbCPU = 1, MaxRAM = 0.25, 
          Kmeans_Only = FALSE, SelectedPCs = FALSE, SpectralFilter = NULL) 
{
  if (!Input_Mask_File == FALSE) {
    driverMask <- get_gdal_info(Input_Mask_File)$driverLongName
    if (driverMask == "ENVI .hdr Labelled") {
      HDR <- read_ENVI_header(get_HDR_name(Input_Mask_File))
      if (!HDR$`data type` == 1) {
        Input_Mask_File <- check_update_mask_format_py(Input_Mask_File, 
                                                    Input_Image_File)
      }
    }
    else {
      Input_Mask_File <- check_update_mask_format_py(Input_Mask_File, 
                                                  Input_Image_File)
    }
  }
  else {
    message("Input_Mask_File not provided in function map_spectral_species.")
    message("Assuming all pixels are valid")
    message("A blank mask will be created for the need of next processing steps")
    HDR <- read_ENVI_header(get_HDR_name(Input_Image_File))
    Mask <- t(matrix(as.integer(1), nrow = HDR$lines, ncol = HDR$samples))
    MaskPath_Update <- paste(file_path_sans_ext(Input_Image_File), 
                             "_BlankMask", sep = "")
    Input_Mask_File <- update_shademask(MaskPath = FALSE, 
                                        HDR = HDR, Mask = Mask, MaskPath_Update = MaskPath_Update)
  }
  Kmeans_info <- NULL
  if (is.null(SpectralSpace_Output$PCA_Files) | is.null(SpectralSpace_Output$TypePCA)) {
    message("Please define input variable SpectralSpace_Output as a list including")
    message("PCA_Files: corresponds to the raster data to be processed (not necessarily resulting from PCA)")
    message("TypePCA: defines main directory where outputs will be written")
    message("This variable is automatically produced as an output of function perform_PCA()")
    message("However, you can set it manually, for example if you want to use spectral indices")
    message("as input raster data instead of PCA file produced from reflectance data")
    stop()
  }
  if (!file.exists(SpectralSpace_Output$PCA_Files)) {
    error_no_PCA_file(SpectralSpace_Output$PCA_Files)
    stop()
  }
  Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, 
                                        SpectralSpace_Output$TypePCA, "SpectralSpecies")
  Output_Dir_PCA <- define_output_subdir(Output_Dir, Input_Image_File, 
                                         SpectralSpace_Output$TypePCA, "PCA")
  Spectral_Species_Path <- file.path(Output_Dir_SS, "SpectralSpecies")
  if (typeof(SelectedPCs) == "logical") {
    if (SelectedPCs == FALSE) {
      PC_Select_Path <- file.path(Output_Dir_PCA, "Selected_Components.txt")
    }
    else {
      message("Error when defining SelectedPCs :")
      message("either set SelectedPCs = FALSE")
      message("or provide a vectorincluding the rank of the variables to be selected from SpectralSpace_Output$PCA_Files")
      stop()
    }
  }
  else {
    PC_Select_Path = "NoFile"
  }
  if (file.exists(PC_Select_Path)) {
    PC_Select <- utils::read.table(PC_Select_Path)[[1]]
  }
  else if (is.numeric(SelectedPCs)) {
    PC_Select <- SelectedPCs
  }
  else {
    error_PC_sel(Output_Dir_PCA)
    stop()
  }
  message("Selected components:")
  print(PC_Select)
  ImNames <- list(Input_Image = Input_Image_File, Mask_list = Input_Mask_File)
  if (is.null(SpectralSpace_Output$nb_partitions)) {
    nb_partitions <- 20
  }
  else {
    nb_partitions <- SpectralSpace_Output$nb_partitions
  }
  Pix_Per_Partition <- define_pixels_per_iter(ImNames, nb_partitions = nb_partitions)
  message("Partition!")
  ImPathHDR <- get_HDR_name(SpectralSpace_Output$PCA_Files)
  HDR <- read_ENVI_header(ImPathHDR)
  message("HDR Read!")
  Subset <- get_random_subset_from_image(ImPath = SpectralSpace_Output$PCA_Files, 
                                         MaskPath = Input_Mask_File, nb_partitions = nb_partitions, 
                                         Pix_Per_Partition = Pix_Per_Partition, kernel = NULL, 
                                         MaxRAM = MaxRAM)
  SubsetInit <- Subset
  dataPCA <- Subset$DataSubset[, PC_Select]
  message("Subset!")
  if (length(PC_Select) == 1) {
    dataPCA <- matrix(dataPCA, ncol = 1)
  }
  print("perform k-means clustering for each subset and define centroids")
  Kmeans_info <- init_kmeans(dataPCA = dataPCA, nb_partitions = nb_partitions, 
                             nbclusters = nbclusters, nbCPU = nbCPU)
  Kmeans_info$SpectralSpecies <- Spectral_Species_Path
  if (Kmeans_info$Error == FALSE) {
    if (Kmeans_Only == FALSE) {
      apply_kmeans(PCA_Path = SpectralSpace_Output$PCA_Files, 
                   PC_Select = PC_Select, Input_Mask_File = Input_Mask_File, 
                   Kmeans_info = Kmeans_info, Spectral_Species_Path = Spectral_Species_Path, 
                   nbCPU = nbCPU, MaxRAM = MaxRAM)
    }
    else {
      print("'Kmeans_Only' was set to TRUE: kmeans was not applied on the full image")
      print("Please set 'Kmeans_Only' to FALSE if you want to produce spectral species map")
    }
    Kmeans_Path <- file.path(Output_Dir_PCA, "Kmeans_Info.RData")
    save(Kmeans_info, file = Kmeans_Path)
  }
  else {
    Output_Dir_Error <- define_output_subdir(Output_Dir, 
                                             Input_Image_File, SpectralSpace_Output$TypePCA, 
                                             "ErrorReport")
    LocError <- unique(c(which(!is.finite(Kmeans_info$MinVal)), 
                         which(!is.finite(Kmeans_info$MaxVal))))
    ValError <- which(!is.finite(dataPCA[, LocError[1]]))
    DataError <- SubsetInit$DataSubset[ValError, ]
    DataErrorCR <- Subset$DataSubset[ValError, ]
    CoordinatesError <- SubsetInit$coordPix[ValError, ]
    FileError <- file.path(Output_Dir_Error, "ErrorPixels.RData")
    ErrorReport <- list(CoordinatesError = CoordinatesError, 
                        DataError = DataError, DataError_afterCR = DataErrorCR, 
                        SpectralFilter = SpectralFilter)
    save(ErrorReport, file = FileError)
    message("")
    message("*********************************************************")
    message("       An error report directory has been produced.      ")
    message("Please check information about data causing errors here: ")
    print(FileError)
    message("               The process will now stop                 ")
    message("*********************************************************")
    message("")
    stop()
  }
  return(Kmeans_info)
}



check_update_mask_format_py <- function (Input_Mask_File, Input_Image_File) 
{
  Mask <- raster(Input_Mask_File)
  ncols <- as.integer(ncol(Mask))
  nrows <- as.integer(nrow(Mask))
  HDR <- read_ENVI_header(get_HDR_name(Input_Image_File))
  if (!ncols == HDR$samples | !nrows == HDR$lines) {
    message("Warning: image and corresponding mask do not have the same dimensions")
    message("Please make sure dimensions match between image and mask")
    stop()
  }
  else {
    MaskPath_Update <- Input_Mask_File
  }
  rm(Mask)
  gc()
  return(MaskPath_Update)
}

diversity_from_plots_nofunc <- function (Raster_SpectralSpecies, Plots, nbclusters = 50, Raster_Functional = FALSE, 
          Selected_Features = FALSE, Name_Plot = FALSE, Hellinger = FALSE, 
          pcelim = 0.02) 
{
  HDR <- read_ENVI_header(paste(Raster_SpectralSpecies, ".hdr", 
                                sep = ""))
  nbRepetitions <- HDR$bands
  nbPlots <- length(Plots)
  Richness.AllRep <- Shannon.AllRep <- Fisher.AllRep <- Simpson.AllRep <- list()
  Richness <- Shannon <- Fisher <- Simpson <- data.frame()
  Raster <- terra::rast(Raster_SpectralSpecies, lyrs = 1)
  XY <- list()
  for (ip in 1:nbPlots) {
    vector_file <- Plots[[ip]]
    if (file.exists(vector_file)) {
      Plot <- terra::vect(vector_file)
      if (!terra::same.crs(Raster, Plot)) 
        Plot <- terra::project(x = Plot, y = Raster)
    }
    else {
      print(paste(vector_file, "cannot be found"))
    }
    XY0 <- extract_pixels_coordinates(x = Raster, y = Plot)
    XY <- c(XY, XY0)
  }
  nbPolygons <- length(XY)
  message(paste("Number of validation plots : ", nbPolygons))
  Pixel_Inventory_All <- Pixel_Hellinger_All <- list()
  pb <- progress::progress_bar$new(format = "computing alpha diversity [:bar] :percent in :elapsed", 
                         total = nbPolygons, clear = FALSE, width = 100)
  for (ip in 1:nbPolygons) {
    pb$tick()
    if (all(is.na(XY[[ip]]$col))) {
      if (length(Name_Plot) == nbPolygons) {
        message(paste("Polygon named", Name_Plot[ip], 
                      "is out of the raster"))
        Name_Plot[ip] <- NA
      }
      Richness <- rbind(Richness, NA, row.names = NULL, 
                        col.names = NULL)
      Fisher <- rbind(Fisher, NA, row.names = NULL, col.names = NULL)
      Shannon <- rbind(Shannon, NA, row.names = NULL, 
                       col.names = NULL)
      Simpson <- rbind(Simpson, NA, row.names = NULL, 
                       col.names = NULL)
    }
    else {
      ExtractIm <- extract.big_raster(Raster_SpectralSpecies, 
                                      XY[[ip]])
      if (length(XY[[ip]]$col) == 1) {
        ExtractIm <- matrix(ExtractIm, ncol = nbRepetitions)
      }
      Pixel_Inventory <- Pixel_Hellinger <- list()
      Richness.tmp <- Shannon.tmp <- Fisher.tmp <- Simpson.tmp <- vector(length = nbRepetitions)
      for (i in 1:nbRepetitions) {
        if (nbRepetitions == 1) {
          Distritab <- table(ExtractIm)
        }
        else {
          Distritab <- table(ExtractIm[, i])
        }
        Pixel_Inventory[[i]] <- as.data.frame(Distritab)
        if (length(which(Pixel_Inventory[[i]]$Var1 == 
                         0)) == 1) {
          Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-which(Pixel_Inventory[[i]]$Var1 == 
                                                                0), ]
        }
        SumPix <- sum(Pixel_Inventory[[i]]$Freq)
        if (SumPix < 25) {
          message("Less than 25 pixels for validation plot")
          if (length(Name_Plot) == nbPolygons) {
            message(Name_Plot[ip])
          }
          message("Please consider applying a buffer")
          message("We recommend at least 25 pixels per plot to compute diversity metrics")
        }
        ThreshElim <- pcelim * SumPix
        ElimZeros <- which(Pixel_Inventory[[i]]$Freq < 
                             ThreshElim)
        if (length(ElimZeros) >= 1) {
          Pixel_Inventory[[i]] <- Pixel_Inventory[[i]][-ElimZeros, 
          ]
        }
        Alpha <- get_alpha_metrics(Pixel_Inventory[[i]]$Freq)
        Richness.tmp[i] <- as.numeric(Alpha$Richness)
        Fisher.tmp[i] <- Alpha$fisher
        Shannon.tmp[i] <- Alpha$Shannon
        Simpson.tmp[i] <- Alpha$Simpson
      }
      Richness.AllRep[[ip]] <- Richness.tmp
      Shannon.AllRep[[ip]] <- Shannon.tmp
      Fisher.AllRep[[ip]] <- Fisher.tmp
      Simpson.AllRep[[ip]] <- Simpson.tmp
      Richness <- rbind(Richness, mean(Richness.tmp), 
                        row.names = NULL, col.names = NULL)
      Fisher <- rbind(Fisher, mean(Fisher.tmp), row.names = NULL, 
                      col.names = NULL)
      Shannon <- rbind(Shannon, mean(Shannon.tmp), row.names = NULL, 
                       col.names = NULL)
      Simpson <- rbind(Simpson, mean(Simpson.tmp), row.names = NULL, 
                       col.names = NULL)
      Pixel_Inventory_All[[ip]] <- Pixel_Inventory
      if (Hellinger == TRUE) {
        for (i in 1:nbRepetitions) {
          Pixel_Hellinger[[i]] <- Pixel_Inventory[[i]]
          Pixel_Hellinger[[i]]$Freq <- sqrt(Pixel_Hellinger[[i]]$Freq/sum(Pixel_Hellinger[[i]]$Freq))
        }
        Pixel_Hellinger_All[[ip]] <- Pixel_Hellinger
      }
    }
  }
  Richness.AllRep <- do.call(rbind, Richness.AllRep)
  Shannon.AllRep <- do.call(rbind, Shannon.AllRep)
  Fisher.AllRep <- do.call(rbind, Fisher.AllRep)
  Simpson.AllRep <- do.call(rbind, Simpson.AllRep)
 
  BC <- list()
  pb <- progress_bar$new(format = "computing beta diversity [:bar] :percent in :elapsed", 
                         total = nbRepetitions, clear = FALSE, width = 100)
  for (i in 1:nbRepetitions) {
    pb$tick()
    MergeDiversity <- matrix(0, nrow = nbclusters, ncol = nbPolygons)
    for (j in 1:nbPolygons) {
      if (nbRepetitions > 1) {
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Inventory_All[[j]][[i]]$Var1))
        SelFrequency <- Pixel_Inventory_All[[j]][[i]]$Freq
      }
      else {
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Inventory_All[[j]][[i]]$ExtractIm))
        SelFrequency <- Pixel_Inventory_All[[j]][[i]]$Freq
      }
      MergeDiversity[SelSpectralSpecies, j] = SelFrequency
    }
    BC[[i]] <- vegan::vegdist(t(MergeDiversity), method = "bray")
  }
  BC_mean <- 0 * BC[[1]]
  for (i in 1:nbRepetitions) {
    BC_mean <- BC_mean + BC[[i]]
  }
  Hellinger_mean <- Hellmat <- NULL
  if (Hellinger == TRUE) {
    Hellmat <- list()
    for (i in 1:nbRepetitions) {
      MergeDiversity <- matrix(0, nrow = nbclusters, ncol = nbPolygons)
      for (j in 1:nbPolygons) {
        SelSpectralSpecies <- as.numeric(as.vector(Pixel_Hellinger_All[[j]][[i]]$Var1))
        SelFrequency <- Pixel_Hellinger_All[[j]][[i]]$Freq
        MergeDiversity[SelSpectralSpecies, j] = SelFrequency
      }
      Hellmat[[i]] <- vegan::vegdist(t(MergeDiversity), 
                                     method = "euclidean")
    }
    Hellinger_mean <- 0 * Hellmat[[1]]
    for (i in 1:nbRepetitions) {
      Hellinger_mean <- Hellinger_mean + Hellmat[[i]]
    }
    Hellinger_mean <- Hellinger_mean/nbRepetitions
  }
  BC_mean <- as.matrix(BC_mean/nbRepetitions)
  names(Richness) <- "Richness"
  names(Fisher) <- "Fisher"
  names(Shannon) <- "Shannon"
  names(Simpson) <- "Simpson"
  return(list(Richness = Richness, Fisher = Fisher, Shannon = Shannon, 
              Simpson = Simpson, fisher.All = Fisher.AllRep, Shannon.All = Shannon.AllRep, 
              Simpson.All = Simpson.AllRep, 
              BCdiss = BC_mean, BCdiss.All = BC, Hellinger = Hellinger_mean, 
              Hellinger_ALL = Hellmat, Name_Plot = Name_Plot))
}

