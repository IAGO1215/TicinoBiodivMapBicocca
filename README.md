# TicinoBiodivMapBicocca

These notebooks and source code are dedicated to calculate biodiversity indices of Ticino Park in Italy, which is the first step of this project. The satellite data are from Sentinel-2 and later from PRISMA. The programming languages are mainly python (jupyter notebooks) and partly R (to utilize the R packages "[preprocS2](https://github.com/jbferet/preprocs2)" and "[biodivMapR](https://github.com/jbferet/biodivMapR)"). 

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
4. Rasterio (and osgeo, which is currently mandatory in order to successfully import the rasterio package.)
5. Matplotlib

## Required R Packages

1. tools
2. raster
3. preprocS2
4. biodivMapR

## Folder structure: 

```
root
├── FieldData
│   ├── Sintesi_siti_ParcoTicino_con_specie.xlsx
│   ├── fieldplots.shp
│   ├── Shapefile
│   │   ├──confini_foreste_finale.shp
│   ├── ShapefileCorretti
│   │   ├──confini_foreste_corretti.shp
│   ├── ShapefileCorretti
│   │   ├──confini_prova.shp
│   ├── Scheda di Campo
├── ProcessedData
│   ├── 2022 06
│   │   ├──2022 06.hdr
│   │   ├──2022 06Cropped.hdr
│   │   ├──2022 06Mask.hdr
│   ├── 2022 09
├── Python
├── R
├── RawData
│   ├── 2022 06
│   ├── 2022 09
├── Results
│   ├── 2022 06Cropped
│   ├── 2022 09Cropped
```

