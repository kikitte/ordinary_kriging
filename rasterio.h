#include "gdal.h"

#ifndef KRIGIN_RASTERIO
#define KRIGIN_RASTERIO

struct RasterInfo
{
  double TOPlEFT_X;
  double TOPLEFT_Y;
  int COLS;
  int ROWS;
  double RESOLUTION;
  double NODATA_VALUE;
  char *GDAL_DRIVER;     // gdal driver for output
  GDALDataType GDAL_GDT; // gdal data type for output
};

CPLErr save_raster(void *raster_array, struct RasterInfo *rast_param, const char *output_filename);

#endif