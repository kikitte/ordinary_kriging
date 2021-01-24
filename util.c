#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>
#include "util.h"

char *read_file_byte(const char *path)
{
  FILE *f = fopen(path, "r");

  if (NULL == f)
  {
    fprintf(stderr, "unable to read file: '%s'\n", path);
    return NULL;
  }

  fseek(f, 0L, SEEK_END);
  long fileSize = ftell(f);

  if (fileSize == 0)
  {
    fclose(f);
    return NULL;
  }

  fseek(f, 0L, SEEK_SET);
  char *fileContent = malloc(sizeof(char) * fileSize + 1);
  char *fileContentPtr = fileContent;
  int ch;
  while ((ch = fgetc(f)) != EOF)
  {
    *fileContentPtr++ = ch;
  }
  *fileContentPtr = '\0';

  fclose(f);
  return fileContent;
}

CPLErr save_raster(void *raster_array, struct RasterInfo *rast_param, const char *output_filename)

{
  GDALAllRegister();

  GDALDriverH hTiff = GDALGetDriverByName(rast_param->GDAL_DRIVER);

  if (NULL == hTiff)
  {
    fprintf(stderr, "unknow driver type: %s\n", rast_param->GDAL_DRIVER);
    return CE_Failure;
  }

  double geoTransform[] = {rast_param->TOPLEFT_Y, rast_param->RESOLUTION, 0, rast_param->TOPlEFT_X, 0, -rast_param->RESOLUTION};

  GDALDatasetH ds = GDALCreate(hTiff, output_filename, rast_param->COLS, rast_param->ROWS, 1, rast_param->GDAL_GDT, NULL);

  GDALSetGeoTransform(ds, geoTransform);

  GDALRasterBandH band1 = GDALGetRasterBand(ds, 1);
  const int COLS = rast_param->COLS, ROWS = rast_param->ROWS;
  const GDALDataType GDT = rast_param->GDAL_GDT;
  CPLErr err = GDALRasterIO(band1, GF_Write, 0, 0, COLS, ROWS, raster_array, COLS, ROWS, GDT, 0, 0);

  GDALClose(ds);

  return err;
}