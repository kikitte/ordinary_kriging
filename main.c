#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"
#include "config.h"
#include "util.h"
#include "cJSON/cJSON.h"

// 根据临近点搜索条件生成sectors
struct SectorsWrap *make_sectorswrap(struct NeighborhoodOption *neighborOpt);
// 分割栅格使其每一部分在不同的线程中运行、计算
struct RasterInfosWrap *split_raster(struct RasterInfo *info);

int main(int argc, char **argv)
{
  if (2 != argc)
  {
    puts("Usage: oridinary_kriging <config>");
    exit(EXIT_FAILURE);
  }
  char *configPath = argv[1];
  cJSON *config = parse_config(configPath);

  char *inputPointsPath = get_input_points_path(config);
  char *outputRasterPath = get_output_raster_path(config);
  struct NeighborhoodOption *neighborOption = get_neighborhood(config);
  struct VariogramModel *variogramOpt = get_variogram(config);
  struct RasterInfo *rasterInfo = get_raster(config);
  struct Points *points = read_input_points(inputPointsPath);

  if (NULL == inputPointsPath || NULL == outputRasterPath ||
      NULL == neighborOption || NULL == variogramOpt ||
      NULL == rasterInfo || NULL == points)
  {
    free(neighborOption);
    free(variogramOpt);
    free(rasterInfo);
    if (points)
    {
      free(points->data);
    }
    free(points);
    cJSON_Delete(config);
    exit(EXIT_FAILURE);
  }

  VariogramFunction modelFunction = getFunctionForVariogramType(variogramOpt->TYPE);
  struct SectorsWrap *sectorsWrap = make_sectorswrap(neighborOption);
  struct RasterInfosWrap *rasterInfosWrap = split_raster(rasterInfo);
  struct OrdinaryKrigingParams ordinaryKrigingParams[rasterInfosWrap->numbers];

  pthread_t threadsId[rasterInfosWrap->numbers];
  for (int i = 0; i < rasterInfosWrap->numbers; ++i)
  {
    ordinaryKrigingParams[i].modelFunction = modelFunction;
    ordinaryKrigingParams[i].points = points;
    ordinaryKrigingParams[i].rasterInfo = rasterInfosWrap->infos + i;
    ordinaryKrigingParams[i].sectorsWrap = sectorsWrap;
    ordinaryKrigingParams[i].variogramOpt = variogramOpt;

    pthread_create(threadsId + i, NULL, ordinary_kriging, ordinaryKrigingParams + i);
  }
  for (int i = 0; i < rasterInfosWrap->numbers; ++i)
  {
    pthread_join(threadsId[i], NULL);
  }

  // 写出栅格
  save_raster(rasterInfo, rasterInfosWrap, outputRasterPath);

  cJSON_Delete(config);
  free(points->data);
  free(points);
  free(sectorsWrap->sectors);
  free(sectorsWrap);
  free(neighborOption);
  free(variogramOpt);
  free(rasterInfo);
  for (int i = 0; i < rasterInfosWrap->numbers; ++i)
  {
    free(rasterInfosWrap->infos[i].rasterArray);
  }
  free(rasterInfosWrap->infos);
  free(rasterInfosWrap);

  return 0;
}

struct SectorsWrap *make_sectorswrap(struct NeighborhoodOption *neighborOpt)
{
  int sectorsCount = 0;
  double sectorInitOffsetAngle = 0;

  if (!strcmp(SECTOR_TYPE_1, neighborOpt->SECTOR_TYPE))
  {
    sectorsCount = 1;
  }
  else if (!strcmp(SECTOR_TYPE_4S, neighborOpt->SECTOR_TYPE))
  {
    sectorsCount = 4;
  }
  else if (!strcmp(SECTOR_TYPE_4S_45D, neighborOpt->SECTOR_TYPE))
  {
    sectorsCount = 4;
    sectorInitOffsetAngle = M_PI_4;
  }
  else if (!strcmp(SECTOR_TYPE_8, neighborOpt->SECTOR_TYPE))
  {
    sectorsCount = 8;
  }
  else
  {
    return NULL;
  }

  struct Sector *sectors = malloc(sizeof(struct Sector) * sectorsCount);

  const double ANGLE_PER_SECTOR = 2 * M_PI / sectorsCount;

  for (int i = 0; i < sectorsCount; ++i)
  {
    struct Sector *sector = &sectors[i];
    sector->minNeighbords = neighborOpt->MIN_NEIGHBORDS;
    sector->maxNeighbords = neighborOpt->MAX_HEIGHBORDS;
    sector->maxDistance = neighborOpt->MAX_DISTANCE;

    if (i)
    {
      struct Sector *lastSector = &sectors[i - 1];
      sector->angleFrom = lastSector->angleTo;
      sector->angleTo = lastSector->angleTo + ANGLE_PER_SECTOR;
    }
    else
    {
      sector->angleFrom = sectorInitOffsetAngle;
      sector->angleTo = sectorInitOffsetAngle + ANGLE_PER_SECTOR;
    }
  }

  struct SectorsWrap *sectorsWrap = malloc(sizeof(struct SectorsWrap));
  sectorsWrap->sectors = sectors;
  sectorsWrap->numbers = sectorsCount;

  return sectorsWrap;
}

struct RasterInfosWrap *split_raster(struct RasterInfo *info)
{
  const int CPU_NUMBERS = get_nprocs();
  const int ROWS = info->ROWS, COLS = info->COLS;
  const int MAX_ROWS_COLS = MAX(ROWS, COLS);

  struct RasterInfosWrap *wrap = malloc(sizeof(struct RasterInfosWrap));

  if (MAX_ROWS_COLS > CPU_NUMBERS)
  {
    int LINES_PER_CPU = MAX_ROWS_COLS / CPU_NUMBERS;

    wrap->numbers = CPU_NUMBERS;
    wrap->infos = malloc(sizeof(struct RasterInfo) * CPU_NUMBERS);

    for (int i = 0; i < CPU_NUMBERS; ++i)
    {
      if (i < CPU_NUMBERS - 1)
      {
        if (ROWS > COLS)
        {
          wrap->infos[i].TOPLEFT_X = info->TOPLEFT_X;
          wrap->infos[i].TOPLEFT_Y = info->TOPLEFT_Y - (info->RESOLUTION * i * LINES_PER_CPU);
          wrap->infos[i].COLS = info->COLS;
          wrap->infos[i].ROWS = LINES_PER_CPU;
        }
        else
        {
          wrap->infos[i].TOPLEFT_X = info->TOPLEFT_X + (info->RESOLUTION * i * LINES_PER_CPU);
          wrap->infos[i].TOPLEFT_Y = info->TOPLEFT_Y;
          wrap->infos[i].COLS = LINES_PER_CPU;
          wrap->infos[i].ROWS = info->ROWS;
        }
      }
      else
      {
        if (ROWS > COLS)
        {
          wrap->infos[i].TOPLEFT_X = info->TOPLEFT_X;
          wrap->infos[i].TOPLEFT_Y = info->TOPLEFT_Y - (info->RESOLUTION * i * LINES_PER_CPU);
          wrap->infos[i].COLS = info->COLS;
          wrap->infos[i].ROWS = ROWS - (i * LINES_PER_CPU);
        }
        else
        {
          wrap->infos[i].TOPLEFT_X = info->TOPLEFT_X + (info->RESOLUTION * i * LINES_PER_CPU);
          wrap->infos[i].TOPLEFT_Y = info->TOPLEFT_Y;
          wrap->infos[i].COLS = COLS - (i * LINES_PER_CPU);
          wrap->infos[i].ROWS = info->ROWS;
        }
      }
    }
  }
  else
  {
    wrap->numbers = MAX_ROWS_COLS;
    wrap->infos = malloc(sizeof(struct RasterInfo) * MAX_ROWS_COLS);

    for (int i = 0; i < MAX_ROWS_COLS; ++i)
    {
      if (ROWS > COLS)
      {
        wrap->infos[i].TOPLEFT_Y = info->TOPLEFT_Y - info->RESOLUTION * i;
        wrap->infos[i].TOPLEFT_X = info->TOPLEFT_X;
        wrap->infos[i].COLS = info->COLS;
        wrap->infos[i].ROWS = 1;
      }
      else
      {
        wrap->infos[i].TOPLEFT_Y = info->TOPLEFT_Y;
        wrap->infos[i].TOPLEFT_X = info->TOPLEFT_X + info->RESOLUTION * i;
        wrap->infos[i].COLS = 1;
        wrap->infos[i].ROWS = info->ROWS;
      }
    }
  }

  for (int i = 0; i < wrap->numbers; ++i)
  {
    wrap->infos[i].RESOLUTION = info->RESOLUTION;
    wrap->infos[i].NODATA_VALUE = info->NODATA_VALUE;
    wrap->infos[i].GDAL_DRIVER = info->GDAL_DRIVER;
    wrap->infos[i].GDAL_GDT = info->GDAL_GDT;
  }

  return wrap;
}