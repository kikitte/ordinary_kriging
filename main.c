#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"
#include "rasterio.h"
#include "config.h"
#include "cJSON/cJSON.h"

// 根据临近点搜索条件生成sectors
struct SectorsWrap *make_sectorswrap(struct NeighborhoodOption *neighborOpt);


struct Points *readSamplePointFromFile(const char *filePath);

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

  if (NULL == inputPointsPath || NULL == outputRasterPath ||
      NULL == neighborOption || NULL == variogramOpt ||
      NULL == rasterInfo)
  {
    free(inputPointsPath);
    free(outputRasterPath);
    free(neighborOption);
    free(variogramOpt);
    free(rasterInfo);
    cJSON_Delete(config);
    exit(EXIT_FAILURE);
  }

  // param1: sample data [x1, y1, v1, x2, y2, v2, ...]
  struct Points *points = readSamplePointFromFile(inputPointsPath);
  if (NULL == points)
  {
    printf("error on reading data from file: %s\n", inputPointsPath);
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

struct Points *readSamplePointFromFile(const char *filePath)
{
  int ch;

  FILE *f = fopen(filePath, "r");
  if (NULL == f)
  {
    return NULL;
  }

  // read the file content
  fseek(f, 0L, SEEK_END);
  long fileByteSize = ftell(f);
  char *fileContent = malloc(fileByteSize + 1);
  char *editableFileContent = fileContent;
  fseek(f, 0L, SEEK_SET);

  while ((ch = getc(f)) != EOF)
  {
    *editableFileContent++ = ch;
  }
  *editableFileContent++ = '\0';

  unsigned lineCount = 0;
  unsigned lastCharIsRN = 0;
  editableFileContent = fileContent;
  while ('\0' != (ch = *editableFileContent++))
  {
    if ('\r' == ch || '\n' == ch)
    {
      *(editableFileContent - 1) = '\n';

      if (!lastCharIsRN)
      {
        ++lineCount;
      }
      lastCharIsRN = 1;
    }
    else
    {
      lastCharIsRN = 0;
    }
  }

  ;
  char *doubleStrEnd;
  unsigned skipThisLine = 0;
  unsigned validPointValueCount = 0;
  double *points = malloc(sizeof(double) * lineCount * 3);
  double *editablePoints = points;

  for (editableFileContent = fileContent, ch = *fileContent; '\0' != ch; ch = *editableFileContent)
  {
    if (skipThisLine)
    {
      if ('\n' == ch)
      {
        skipThisLine = 0;
      }
      ++editableFileContent;
      continue;
    }

    if ('\n' == ch && '\n' == *(editableFileContent + 1))
    {
      ++editableFileContent;
      continue;
    }

    double val = strtod(editableFileContent, &doubleStrEnd);
    if (editableFileContent == doubleStrEnd)
    {
      skipThisLine = 1;
      ++editableFileContent;
      continue;
    }
    else
    {
      ++validPointValueCount;
      *editablePoints++ = val;
      editableFileContent = doubleStrEnd + 1;
    }
  }

  double *validPoints = malloc(sizeof(double) * validPointValueCount);
  double *ediatbleValidPoints = validPoints;

  editablePoints = points;
  for (int i = 0; i < validPointValueCount; ++i)
  {
    *ediatbleValidPoints++ = *editablePoints++;
  }

  free(points);
  free(fileContent);
  fclose(f);

  // *pointNumbers = validPointValueCount / 3;
  struct Points *pointsStrucure = malloc(sizeof(struct Points));
  pointsStrucure->numbers = validPointValueCount / 3;
  pointsStrucure->data = validPoints;
  return pointsStrucure;
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
