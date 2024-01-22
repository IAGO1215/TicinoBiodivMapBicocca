# TicinoBiodivMapBicocca

These notebooks and source code are dedicated to calculate biodiversity indices of Ticino Park in Italy based on the satellite data (2022), which is the first step of this project. The satellite data are from Sentinel-2 and later from PRISMA. The programming languages are mainly python (jupyter notebooks) and partly R (to utilize the R packages "[preprocS2](https://github.com/jbferet/preprocs2)" and "[biodivMapR](https://github.com/jbferet/biodivMapR)"). 

The next and also last step is to validate the calculation results based on the ground truth data collected from field survey (2022). 

## Workflow

1. Obtain raw satellite data and save them in different subfolders separate by dates (months).
2. Process these raw satellite data in R Studio and create corresponding stacked raster files.
3. Crop the stacked raster files to the minimum extent which covers the whole RoI.
4. Process the field shapefile and create mask layer(s) from it.
5. Use the cropped stacked raster files and mask layers into R Studio and utilize "biodivMapR" package to calculate biodiversity indices.
6. Calculate biodiversity indices based on field data (ground truth).
7. Perform validation. 

## Required Python Packages

1. Numpy
2. Pandas
3. Shapely
4. GeoPandas
5. Rasterio (and osgeo, which is currently mandatory in order that the "rasterio" package can be successfully imported.)
6. Matplotlib
7. Seaborn

## Required R Packages

1. tools
2. raster
3. preprocS2
4. biodivMapR

## Folder structure: 

```
root
├── FieldData
│   ├── Shapefile (THE ORIGINAL SHAPEFILE OF TICINO PARK)
│   │   ├──confini_foreste_finale.shp
│   ├── ShapefileAggiunto (THE ADDITIONAL SHAPEFILE OF TICINO PARK WHICH ADDS SOME FOREST AREAS)
│   │   ├──ticino_classification.shp
│   ├── ShapefileCorretti (THE ORIGINAL SHAPEFILE OF TICINO PARK CORRECTED IN PYTHON)
│   │   ├──confini_foreste_corretti.shp
│   ├── ShapefileCorrettiAggiunto (THE FINAL MERGED CORRECTED SHAPEFILE OF TICINO PARK INCLUDING ALL THE ROI)
│   │   ├──confini_foreste_estesi.shp
├── PRISMA Classification Indices (THE STACKED VEGETATION INDICES RASTERS CREATED FROM PRISMA RASTERS)
├── PRISMA Raster Mask (THE MASK RASTERS FOR PRISMA RASTERS)
├── PRISMA Raster Raw (THE ORIGINAL PRISMA RASTERS)
│   ├── Merged (MERGED PRISMA RASTERS)
│   ├── Original (RAW PRISMA RASTERS)
├── ProcessedData (SENTINEL-2 STACKED BAND DATA)
│   ├── 2022 06
│   │   ├──2022 06.hdr
│   │   ├──2022 06Cropped.hdr
│   │   ├──2022 06Mask.hdr
│   ├── 2022 09
├── Python
├── R
├── RawData (SENTINEL-2 RAW DATA)
│   ├── 2022 06
│   ├── 2022 09
├── Results (SENTINEL-2 BIODIVMAPR RESULTS)
│   ├── 2022 06Cropped
│   ├── 2022 09Cropped
├── ResultsPRISMA (PRISMA BIODIVMAPR RESULTS)
```

