#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "ordinary_kriging.h"
#include "variogram.h"

// #define DEBUG
// #define BINNING_DEBUG

double *readSamplePointFromFile(char *, int *);

// 区间分组：用于计算适配模型的参数
void binning(struct Points *points, double lag, int lagNumbers, int *lagsSemivarCount, double *lagsSemivarSum);
// 在没有制定输出栅格参数时根据已有条件生成默认参数
struct RasterInfo *deafultRasterFromPoints(struct Points *points);

int main(int argc, char **argv)
{
  // process args
  if (9 != argc)
  {
    puts("missing arguments.");
    puts("kriging <sample_points> <variogram_model> <lag> <lag_numbers> <min_neighbors> <max_neighbors> <sector_type> <raster_cell_size>");
    exit(EXIT_FAILURE);
  }
  char *argSamplePoints = argv[1];
  char *argVariogramModel = argv[2];
  char *argLag = argv[3];
  char *argLagNumbers = argv[4];
  char *argMinNeighbors = argv[5];
  char *argMaxNeighbors = argv[6];
  char *argSectorType = argv[7];
  char *argRasterCellSize = argv[8];

  // param1: sample data [x1, y1, v1, x2, y2, v2, ...]
  int samplePointNumbers = 0;
  double *samplePoints = readSamplePointFromFile(argSamplePoints, &samplePointNumbers);
  if (NULL == samplePoints)
  {
    printf("error on reading data from file: %s\n", argSamplePoints);
    exit(EXIT_FAILURE);
  }

  // param2: variogram model
  const int variogramModel = atoi(argVariogramModel);
  if (variogramModel != VARIOGRAM_MODEL_LINEAR &&
      variogramModel != VARIOGRAM_MODEL_SPHERICAL &&
      variogramModel != VARIOGRAM_MODEL_EXPONENTIAL &&
      variogramModel != VARIOGRAM_MODEL_GAUSSIAN)
  {
    printf("unknow variogram model: %s\n", argVariogramModel);
    exit(EXIT_FAILURE);
  }

  // param3 & 4: lag & lag numbers
  const double lag = atof(argLag);
  const int lagNumbers = atoi(argLagNumbers);
  if (lag <= 0)
  {
    printf("invalid lag: %s\n", argLag);
    exit(EXIT_FAILURE);
  }
  if (lagNumbers <= 0)
  {
    printf("invalid lag numbers: %s\n", argLagNumbers);
    exit(EXIT_FAILURE);
  }

  // param5 & 6: min & max neighborhoods
  const int minNeighbors = atoi(argMinNeighbors);
  const int maxNeighbors = atoi(argMaxNeighbors);
  if (minNeighbors <= 0)
  {
    printf("invalid minimum neighbords: %s\n", argMinNeighbors);
    exit(EXIT_FAILURE);
  }
  if (maxNeighbors <= 0)
  {
    printf("invalid maximum neighbords: %s\n", argMaxNeighbors);
    exit(EXIT_FAILURE);
  }

  // param 7: sector type
  const int sectorType = atoi(argSectorType);
  if (sectorType <= 0)
  {
    printf("invalid sector type: %s\n", argSectorType);
    exit(EXIT_FAILURE);
  }

  // param8: raster cell size
  const double rasterCellSize = atof(argRasterCellSize);
  if (rasterCellSize <= 0)
  {
    printf("invalid lag raster cell size: %s\n", argRasterCellSize);
    exit(EXIT_FAILURE);
  }

  // ordinaryKriging(samplePoints, samplePointNumbers, variogramModel, lag, lagNumbers, minNeighbors, maxNeighbors, sectorType, rasterCellSize);
  struct Points krigingPoints = {samplePoints, samplePointNumbers};
  struct VariogramModel variogramOpt = {.VAR = VARIOGRAM_MODEL_SPHERICAL,
                                        .C0 = 7.179056886219061,
                                        .CX = 7.179056886219061 + 70.42557944918185,
                                        .A = 480000};
  struct NeighborhoodOption neighborOption = {sectorType, 2, 5, 480000};
  struct RasterInfo rasterInfo = {
      .top = 881249.623773591,
      .left = 227892.999058652,
      .resolution = 200,
      .cols = 2694,
      .rows = 4064,
      .nodata = 0};

  // 2. 根据需要生成半变异图参数：nugget(ArcGIS使用weighted least squares算法生成), sill, range等
  //    根据需要生成lag & lagNumbers

  // int lagsSemivarCount[pointNumbers];
  // double lagsSemivarSum[pointNumbers];
  // struct PointPairDistance pointPairDistance[(int)(pointNumbers * (pointNumbers - 1) / 2)]; // 点与点之间配对个数
  // binning(points, semivarOpt, lagsSemivarCount, lagsSemivarSum, pointPairDistance);

  ordinaryKriging(&krigingPoints, &variogramOpt, &neighborOption, &rasterInfo);

  return 0;
}

double *readSamplePointFromFile(char *filePath, int *pointNumbers)
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
  char *fileContent = malloc(fileByteSize);
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

#ifdef DEBUG
  for (int i = 0, len = validPointValueCount / 3; i < len; ++i)
  {
    printf("%f, %f, %f\n", validPoints[i * 3], validPoints[i * 3 + 1], validPoints[i * 3 + 2]);
  }
#endif

  free(points);
  free(fileContent);
  fclose(f);

  *pointNumbers = validPointValueCount / 3;
  return validPoints;
}

void binning(struct Points *points,
             const double lag,
             const int lagNumbers,
             int *lagsSemivarCount,
             double *lagsSemivarSum)
{
  const double *pointsData = points->data;
  const int pointNumbers = points->numbers;

  for (int i = 0; i < lagNumbers; ++i)
  {
    // 初始化
    lagsSemivarCount[i] = 0;
    lagsSemivarSum[i] = 0;
  }
  for (int i = 0; i < pointNumbers; ++i)
  {
    const int ipIndex = i * 3;
    double ipx = pointsData[ipIndex];
    double ipy = pointsData[ipIndex + 1];
    double ipv = pointsData[ipIndex + 2];
    for (int j = i + 1; j < pointNumbers; ++j)
    {
      const int jpIndex = j * 3;
      double jpx = pointsData[jpIndex];
      double jpy = pointsData[jpIndex + 1];
      double jpv = pointsData[jpIndex + 2];

      const double ijDistance = sqrt(pow(ipx - jpx, 2) + pow(ipy - jpy, 2));
      // 如果lagsSemivarIndex小于lagNumbers则表明距离在lag * lagNumbers以内。
      const int lagsSemivarIndex = ijDistance / lag;
      if (lagsSemivarIndex < lagNumbers)
      {
        const double semivariance = pow(ipv - jpv, 2) / 2;
        ++lagsSemivarCount[lagsSemivarIndex];
        lagsSemivarSum[lagsSemivarIndex] += semivariance;
      }
    }
  }

#ifdef BINNING_DEBUG
  printf("x = [");
  for (int i = 0; i < lagNumbers; ++i)
  {
    if (lagsSemivarCount[i])
    {
      printf("%d%s", (int)(lag * (i + 0.5)), i < lagNumbers - 1 ? ", " : "");
    }
  }
  printf("]\n");

  printf("y = [");
  for (int i = 0; i < lagNumbers; ++i)
  {
    if (lagsSemivarCount[i])
    {
      printf("%f%s", lagsSemivarSum[i] / lagsSemivarCount[i], i < lagNumbers - 1 ? ", " : "");
    }
  }
  printf("]\n");

#endif
}


struct RasterInfo *deafultRasterFromPoints(struct Points *points)
{
  struct RasterInfo *info = malloc(sizeof(struct RasterInfo));

  double minx = DBL_MAX, miny = DBL_MAX, maxx = DBL_MIN, maxy = DBL_MIN;
  double *pointsData = points->data;
  const int pointNumbers = points->numbers;
  for (int i = 0; i < pointNumbers; ++i)
  {
    int ptIndex = 3 * i;
    double x = pointsData[ptIndex];
    double y = pointsData[ptIndex + 1];

    if (x < minx)
    {
      minx = x;
    }
    if (x > maxx)
    {
      maxx = x;
    }
    if (y < miny)
    {
      miny = y;
    }
    if (y > maxy)
    {
      maxy = y;
    }
  }

  info->top = maxy;
  info->left = minx;
  // 默认x方向像元数量为250个
  info->cols = 250;
  info->resolution = (int)((maxx - minx) / info->cols);
  info->rows = (int)((maxy - miny) / info->resolution) + 1;
  return info;
}