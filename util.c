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

CPLErr save_raster(struct RasterInfo *rasterInfo, struct RasterInfosWrap *rasterPartsWrap, const char *output_filename)

{
  GDALAllRegister();

  GDALDriverH hTiff = GDALGetDriverByName(rasterInfo->GDAL_DRIVER);

  if (NULL == hTiff)
  {
    fprintf(stderr, "unknow driver type: %s\n", rasterInfo->GDAL_DRIVER);
    return CE_Failure;
  }

  double geoTransform[] = {rasterInfo->TOPLEFT_X, rasterInfo->RESOLUTION, 0, rasterInfo->TOPLEFT_Y, 0, -rasterInfo->RESOLUTION};

  GDALDatasetH ds = GDALCreate(hTiff, output_filename, rasterInfo->COLS, rasterInfo->ROWS, 1, rasterInfo->GDAL_GDT, NULL);

  GDALSetGeoTransform(ds, geoTransform);

  GDALRasterBandH band1 = GDALGetRasterBand(ds, 1);
  const double RESOLUTION = rasterInfo->RESOLUTION;
  const GDALDataType GDT = rasterInfo->GDAL_GDT;
  for (int i = 0; i < rasterPartsWrap->numbers; ++i)
  {
    struct RasterInfo *partInfo = rasterPartsWrap->infos + i;
    const int OFFSET_COLS = round((partInfo->TOPLEFT_X - rasterInfo->TOPLEFT_X) / RESOLUTION),
              OFFSET_ROWS = round((rasterInfo->TOPLEFT_Y - partInfo->TOPLEFT_Y) / RESOLUTION);
    const int COLS = partInfo->COLS, ROWS = partInfo->ROWS;
    CPLErr err = GDALRasterIO(band1, GF_Write, OFFSET_COLS, OFFSET_ROWS, COLS, ROWS, partInfo->rasterArray, COLS, ROWS, GDT, 0, 0);
  }

  GDALClose(ds);

  return CE_None;
}