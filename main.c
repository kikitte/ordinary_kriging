#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"
#include "config.h"
#include "util.h"
#include "cJSON/cJSON.h"

// 根据临近点搜索条件生成sectors
struct SectorsWrap *make_sectorswrap(struct NeighborhoodOption *neighborOpt);

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
    if (points){
      free(points->data);
    }
    free(points);
    cJSON_Delete(config);
    exit(EXIT_FAILURE);
  }

  VariogramFunction modelFunction = getFunctionForVariogramType(variogramOpt->TYPE);
  struct SectorsWrap *sectorsWrap = make_sectorswrap(neighborOption);

  void *rasterArray = ordinary_kriging(points, variogramOpt, sectorsWrap, rasterInfo, modelFunction);

  // 写出栅格
  save_raster(rasterArray, rasterInfo, outputRasterPath);

  cJSON_Delete(config);
  free(points->data);
  free(points);
  free(sectorsWrap->sectors);
  free(sectorsWrap);
  free(rasterArray);
  free(neighborOption);
  free(variogramOpt);
  free(rasterInfo);

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
  sectorsWrap->count = sectorsCount;

  return sectorsWrap;
}
