#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ordinary_kriging.h"
#include "variogram.h"

#define BINNING_DEBUG

struct RasterInfo *deafultRasterFromPoints(struct Points *points);
// 根据扇区类型搜索用于克里金插值的点
void searchNeighborhoods(struct PointsIndex *neighbords, struct Points *points, struct PointPairDistance *pointPairDistance, struct NeighborhoodOption *neighborOpt);
void binning(struct Points *points, struct SemivariogramOption *semivarOpt, int *lagsSemivarCount, double *lagsSemivarSum, struct PointPairDistance *pointPairDistance);

double *ordinaryKriging(struct Points *points,
                        struct SemivariogramOption *semivarOpt,
                        struct NeighborhoodOption *neighborOpt,
                        struct RasterInfo *rasterInfo)
{
  // 默认生成的RasterInfo结构体通过动态请求内存, 所以后面需要释放
  int isDefaultRasterInfo = 0;
  if (rasterInfo == NULL)
  {
    isDefaultRasterInfo = 1;
    rasterInfo = deafultRasterFromPoints(points);
  }

  // 1. 区间分组，忽略数据的各项异性。
  const int pointNumbers = points->numbers;
  int lagsSemivarCount[pointNumbers];
  double lagsSemivarSum[pointNumbers];
  struct PointPairDistance pointPairDistance[(int)(pointNumbers * (pointNumbers - 1) / 2)]; // 点与点之间配对个数
  binning(points, semivarOpt, lagsSemivarCount, lagsSemivarSum, pointPairDistance);

  // 2. 根据需要生成半变异图参数：nugget(ArcGIS使用weighted least squares算法生成), sill, range等
  //    根据需要生成lag & lagNumbers

  // 3. 栅格点插值
  double *raster = malloc(sizeof(double) * rasterInfo->xsize * rasterInfo->ysize); // 栅格以行存储
  double *rasterPointer = raster;
  for (int x = 0; x < rasterInfo->xsize; ++x)
  {
    double ptx = rasterInfo->left + rasterInfo->resolution * (x + 0.5);
    for (int y = 0; y < rasterInfo->ysize; ++y)
    {
      double pty = rasterInfo->top - rasterInfo->resolution * (y + 0.5);

      double result = 0;

      //

      *rasterPointer++ = result;
    }
  }

  if (isDefaultRasterInfo)
  {
    free(rasterInfo);
  }
  return NULL;
}

void binning(struct Points *points,
             struct SemivariogramOption *semivarOpt,
             int *lagsSemivarCount,
             double *lagsSemivarSum,
             struct PointPairDistance *pointPairDistance)
{
  const double *pointsData = points->data;
  const int pointNumbers = points->numbers;
  const double lag = semivarOpt->lag;
  const int lagNumbers = semivarOpt->lagNumbers;

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
      struct PointPairDistance *pointPair = pointPairDistance++;
      pointPair->fromIndex = i;
      pointPair->toIndex = j;
      pointPair->distance = ijDistance;
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
  info->xsize = 250;
  info->resolution = (int)((maxx - minx) / info->xsize);
  info->ysize = (int)((maxy - miny) / info->resolution) + 1;
  return info;
}

void searchNeighborhoods(struct PointsIndex *neighbords, struct Points *points, struct PointPairDistance *pointPairDistance, struct NeighborhoodOption *neighborOpt)
{
  const double distance = neighborOpt->maxDistance < 0 ? DBL_MAX : neighborOpt->maxDistance;
}