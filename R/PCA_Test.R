###### PCA TEST ######

### Sep NEW ###

### The Whole Image ### 
test_RasterIndice <- file.path(path_abs,"PRISMA Classification Indices","Sel","PRS_L2D_STD_20220906_20220911_NS_mosaic_crop_smooth_v2i_new_StackedIndices_20VL")
test_RasterIndice_MaskOriginal <- file.path(path_abs,"PRISMA Raster Mask","PRS_L2D_STD_20220906_20220911_NS_mosaic_crop_smooth_v2i Mask Original")

temp_RasterIndice_PCA <- biodivMapR::perform_PCA(Input_Image_File = test_RasterIndice,
                        Output_Dir = file.path(outDir_PRISMA,"New Test","Whole"),
                        Input_Mask_File = test_RasterIndice_MaskOriginal,
                        TypePCA = typePCA,
                        FilterPCA = filterPCA,
                        nbCPU = nbCPU,
                        MaxRAM = MaxRAM,
                        Continuum_Removal = FALSE, 
                        Excluded_WL = spectral_excluded)
print(paste0('PCA generated for ',test_RasterIndice))

temp_Summary <- temp_RasterIndice_PCA$PCA_model$sdev
temp_Summary_df <- data.frame(temp_Summary)
names(temp_Summary_df)  = c("Standard Deviation")
write.table(temp_Summary, file = file.path(outDir_PRISMA,"New Test","Whole","PCA R.csv"),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
rm(temp_RasterIndice_PCA,temp_Summary,temp_Summary_df)

### The RoI ###

test_RasterIndice <- file.path(path_abs,"PRISMA Classification Indices","Sel","PRS_L2D_STD_20220906_20220911_NS_mosaic_crop_smooth_v2i_new_StackedIndices_20VL")

temp_RasterIndice_PCA <- biodivMapR::perform_PCA(Input_Image_File = test_RasterIndice,
                        Output_Dir = file.path(outDir_PRISMA,"New Test","RoI"),
                        Input_Mask_File = path_PRISMAMask[[2]],
                        TypePCA = typePCA,
                        FilterPCA = filterPCA,
                        nbCPU = nbCPU,
                        MaxRAM = MaxRAM,
                        Continuum_Removal = FALSE, 
                        Excluded_WL = spectral_excluded)
print(paste0('PCA generated for ',test_RasterIndice))

temp_Summary <- temp_RasterIndice_PCA$PCA_model$sdev
temp_Summary_df <- data.frame(temp_Summary)
names(temp_Summary_df)  = c("Standard Deviation")
write.table(temp_Summary, file = file.path(outDir_PRISMA,"New Test","RoI","PCA R.csv"),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
rm(temp_RasterIndice_PCA,temp_Summary,temp_Summary_df)

### June-July ###

### The Whole Image ### 
test_RasterIndice <- file.path(path_abs,"PRISMA Classification Indices","Sel","PRS_L2D_STD_20220611_20220710_NS_mosaic_crop_smooth_v2i_StackedIndices_20VL")
test_RasterIndice_MaskOriginal <- file.path(path_abs,"PRISMA Raster Mask","PRS_L2D_STD_20220611_20220710_NS_mosaic_crop_smooth_v2i Mask Original")

temp_RasterIndice_PCA <- biodivMapR::perform_PCA(Input_Image_File = test_RasterIndice,
                                                 Output_Dir = file.path(outDir_PRISMA,"New Test","Whole"),
                                                 Input_Mask_File = test_RasterIndice_MaskOriginal,
                                                 TypePCA = typePCA,
                                                 FilterPCA = filterPCA,
                                                 nbCPU = nbCPU,
                                                 MaxRAM = MaxRAM,
                                                 Continuum_Removal = FALSE, 
                                                 Excluded_WL = spectral_excluded)
print(paste0('PCA generated for ',test_RasterIndice))

temp_Summary <- temp_RasterIndice_PCA$PCA_model$sdev
temp_Summary_df <- data.frame(temp_Summary)
names(temp_Summary_df)  = c("Standard Deviation")
write.table(temp_Summary, file = file.path(outDir_PRISMA,"New Test","Whole","PCA R.csv"),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
rm(temp_RasterIndice_PCA,temp_Summary,temp_Summary_df)

### The RoI ###

test_RasterIndice <- file.path(path_abs,"PRISMA Classification Indices","Sel","PRS_L2D_STD_20220611_20220710_NS_mosaic_crop_smooth_v2i_StackedIndices_20VL")

temp_RasterIndice_PCA <- biodivMapR::perform_PCA(Input_Image_File = test_RasterIndice,
                                                 Output_Dir = file.path(outDir_PRISMA,"New Test","RoI"),
                                                 Input_Mask_File = path_PRISMAMask[[1]],
                                                 TypePCA = typePCA,
                                                 FilterPCA = filterPCA,
                                                 nbCPU = nbCPU,
                                                 MaxRAM = MaxRAM,
                                                 Continuum_Removal = FALSE, 
                                                 Excluded_WL = spectral_excluded)
print(paste0('PCA generated for ',test_RasterIndice))

temp_Summary <- temp_RasterIndice_PCA$PCA_model$sdev
temp_Summary_df <- data.frame(temp_Summary)
names(temp_Summary_df)  = c("Standard Deviation")
write.table(temp_Summary, file = file.path(outDir_PRISMA,"New Test","RoI","PCA R.csv"),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)
rm(temp_RasterIndice_PCA,temp_Summary,temp_Summary_df)