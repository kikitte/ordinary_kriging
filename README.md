## Ordinary Kriging Interpolation

This is an c language implemetation using arcgis searching strategy.

### 1. Highlight

- Use arcgis searching strategy to search neighbords of pixel. It produce the same result if semivariogram and it's parameter are the same.
- Thread level paralle. The speed is 4x than ArcGIS.
- Configuring in json file.

### 2. How to use

You might install GDAL in your computer, but as it may be different among diffent system, so we ignore it here.

```bash
# clone soure code
git clone https://www.github.com/kikitte/ordinary_kriging
# compiling
cd ordinary_kriging
make
# run with testing data
build/ordinary_kriging test/config.json
```

Configuration introduction

```
{
  // input: the path of sample points
  "INPUT_POINTS": "/home/kikitte/geoanalyst/kriging/chap15/stations.json",
  // output: the path of result raster
  "OUTPUT_RASTER": "/home/kikitte/geoanalyst/kriging/debug/out2.tif",
  // semivariogram model configuration
  "VariogramModel": {
    // the type of semivarogram model
    "TYPE": "Spherical",
    // nugget
    "NUGGET": 7.179056886219,
    // sill
    "SILL": 77.604636335401,
    // range
    "RANGE": 480000
  },
  // Searching strategy configuration
  "NeighborhoodOption": {
    // The type of sector："OneSector", "FourSector", "FourSector45Offset", "EightSector"
    "SECTOR_TYPE": "FourSector45Offset",
    // The mininum number in each sector
    "MIN_NEIGHBORDS": 2,
    // The maxinum number in each sector
    "MAX_HEIGHBORDS": 5,
    // The maxinum searching distance of each sector
    "MAX_DISTANCE": 480000
  },
  // Output raster configuration
  "RasterInfo": {
    // x coordinate of raster top-left corner
    "TOPLEFT_X": 227942.9990586524945684,
    // y coordinate of raster top-left corner
    "TOPLEFT_Y": 881299.6237735910108313,
    // raster pixel size
    "RESOLUTION": 100,
    // the column number of raster
    "COLS": 5388,
    // the row number of raster
    "ROWS": 8128,
    // nodata value of raster cell
    "NODATA_VALUE": 0,
    // raster output format. reference: https://gdal.org/drivers/raster/index.html
    "GDAL_DRIVER": "GTiff",
    // raster pixel type
    "GDAL_GDT": "Float32"
  }
}
```

### 3. License
### 3. Credit

None

### 4. Thinks

程序的完成参考或使用了很多其他人的代码，在此一并感谢，如下：
This project use a lot of codes from internet, thanks all: 

- [cJSON](https://github.com/DaveGamble/cJSON) 

  for parsing configuration file and reading sample point file.

- solve Ax=B

  [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode)

  [LU decompossion](https://en.wikipedia.org/wiki/LU_decomposition#C_code_example)

  [Don’t invert that matrix](https://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/)

  [Triangular_matrix Forward_and_back_substitution](https://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_back_substitution)

  There is no need to invert the hole matrix and than to solve Ax=B, using gaussian-dissolve and bask-substitution is enough.

- Seimvariogram Model

  This repository code is used: https://github.com/cnca/Kriging

- Esri or ArcGIS
  interesting.