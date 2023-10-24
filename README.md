# TicinoBiodivMapBicocca

These notebooks and source code are dedicated to calculate biodiversity indices of Ticino Park in Italy, which is the first step of this project. The satellite data are from Sentinel-2 and later from PRISMA. The programming languages are mainly python (jupyter notebooks) and partly R (to utilize the R packages "preprocS2" and "biodivMapR"). 

The last step is to validate the calculation results based on the ground truth data collected from field survey. 

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
3. GeoPandas
4. Rasterio (and osgeo)
5. Matplotlib

## Folder structure: 

```

root

```

