#include <stdio.h>
#include <string.h>
#include <gdal.h>
#include "config.h"

cJSON *parse_config(const char *configPath)
{
  FILE *f = fopen(configPath, "r");

  if (NULL == f)
  {
    fprintf(stderr, "unable to read file: '%s'\n", configPath);
    return NULL;
  }

  // read the file content
  fseek(f, 0L, SEEK_END);
  long fileSize = ftell(f);
  fseek(f, 0L, SEEK_SET);
  char *fileContent = malloc(sizeof(char) * fileSize + 1);
  char *fileContentPtr = fileContent;
  int ch;
  while ((ch = fgetc(f)) != EOF)
  {
    *fileContentPtr++ = ch;
  }
  *fileContentPtr++ = '\0';

  cJSON *config = cJSON_Parse(fileContent);
  free(fileContent);

  if (NULL == config)
  {
    const char *error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL)
    {
      fprintf(stderr, "Error before: %s\n", error_ptr);
    }
    return NULL;
  }
  else
  {
    return config;
  }
}

struct VariogramModel *get_variogram(cJSON *config)
{
  cJSON *variogramConf = cJSON_GetObjectItem(config, "VariogramModel");

  if (NULL == variogramConf)
  {
    fprintf(stderr, "missing 'VariogramModel' property\n");
    return NULL;
  }

  struct VariogramModel *vairogramModel = malloc(sizeof(struct VariogramModel));
  vairogramModel->TYPE = cJSON_GetObjectItem(variogramConf, "TYPE")->valuestring;
  vairogramModel->NUGGET = cJSON_GetObjectItem(variogramConf, "NUGGET")->valuedouble;
  vairogramModel->RANGE = cJSON_GetObjectItem(variogramConf, "RANGE")->valuedouble;
  vairogramModel->SILL = cJSON_GetObjectItem(variogramConf, "SILL")->valuedouble;

  return vairogramModel;
}
struct NeighborhoodOption *get_neighborhood(cJSON *config)
{
  cJSON *neighborhoodConf = cJSON_GetObjectItem(config, "NeighborhoodOption");
  if (NULL == neighborhoodConf)
  {
    fprintf(stderr, "missing 'NeighborhoodOption' property\n");
    return NULL;
  }

  struct NeighborhoodOption *neighborhoodOption = malloc(sizeof(struct NeighborhoodOption));
  neighborhoodOption->MAX_DISTANCE = cJSON_GetObjectItem(neighborhoodConf, "MAX_DISTANCE")->valuedouble;
  neighborhoodOption->MAX_HEIGHBORDS = cJSON_GetObjectItem(neighborhoodConf, "MAX_HEIGHBORDS")->valueint;
  neighborhoodOption->MIN_NEIGHBORDS = cJSON_GetObjectItem(neighborhoodConf, "MIN_NEIGHBORDS")->valueint;
  neighborhoodOption->SECTOR_TYPE = cJSON_GetObjectItem(neighborhoodConf, "SECTOR_TYPE")->valuestring;

  return neighborhoodOption;
}
struct RasterInfo *get_raster(cJSON *config)
{
  cJSON *rasterConf = cJSON_GetObjectItem(config, "RasterInfo");
  if (NULL == rasterConf)
  {
    fprintf(stderr, "missing 'RasterInfo' property\n");
    return NULL;
  }

  struct RasterInfo *rasterInfo = malloc(sizeof(struct RasterInfo));
  rasterInfo->TOPlEFT_X = cJSON_GetObjectItem(rasterConf, "TOPlEFT_X")->valuedouble;
  rasterInfo->TOPLEFT_Y = cJSON_GetObjectItem(rasterConf, "TOPLEFT_Y")->valuedouble;
  rasterInfo->RESOLUTION = cJSON_GetObjectItem(rasterConf, "RESOLUTION")->valuedouble;
  rasterInfo->COLS = cJSON_GetObjectItem(rasterConf, "COLS")->valueint;
  rasterInfo->ROWS = cJSON_GetObjectItem(rasterConf, "ROWS")->valueint;
  rasterInfo->NODATA_VALUE = cJSON_GetObjectItem(rasterConf, "NODATA_VALUE")->valuedouble;
  rasterInfo->GDAL_DRIVER = cJSON_GetObjectItem(rasterConf, "GDAL_DRIVER")->valuestring;

  cJSON *gdtItem = cJSON_GetObjectItem(rasterConf, "GDAL_GDT");
  char *gdalDataType;
  if (gdtItem)
  {
    gdalDataType = gdtItem->valuestring;
  }
  else
  {
    fprintf(stderr, "missing 'RasterInfo.GDAL_GDT' property\n");
    free(rasterInfo);
    return NULL;
  };

  if (!strcmp(gdalDataType, "Byte"))
  {
    rasterInfo->GDAL_GDT = GDT_Byte;
  }
  else if (!strcmp(gdalDataType, "UInt16"))
  {
    rasterInfo->GDAL_GDT = GDT_UInt16;
  }
  else if (!strcmp(gdalDataType, "Int16"))
  {
    rasterInfo->GDAL_GDT = GDT_Int16;
  }
  else if (!strcmp(gdalDataType, "UInt32"))
  {
    rasterInfo->GDAL_GDT = GDT_UInt32;
  }
  else if (!strcmp(gdalDataType, "Int32"))
  {
    rasterInfo->GDAL_GDT = GDT_Int32;
  }
  else if (!strcmp(gdalDataType, "Float32"))
  {
    rasterInfo->GDAL_GDT = GDT_Float32;
  }
  else if (!strcmp(gdalDataType, "Float64"))
  {
    rasterInfo->GDAL_GDT = GDT_Float64;
  }
  else
  {
    fprintf(stderr, "Unknown GDAL_GDT");
    free(rasterInfo);
    return NULL;
  }

  return rasterInfo;
}

char *get_input_points_path(cJSON *config)
{
  cJSON *item = cJSON_GetObjectItem(config, "INPUT_POINTS");
  if (item)
  {
    return item->valuestring;
  }
  else
  {
    fprintf(stderr, "missing 'INPUT_POINTS' property\n");
    return NULL;
  }
}

char *get_output_raster_path(cJSON *config)
{
  cJSON *item = cJSON_GetObjectItem(config, "OUTPUT_RASTER");
  if (item)
  {
    return item->valuestring;
  }
  else
  {
    fprintf(stderr, "missing 'OUTPUT_RASTER' property\n");
    return NULL;
  }
}