#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"

// #define BINNING_DEBUG
#define MAKE_SECTORS_DEBUG

// 在没有制定输出栅格参数时根据已有条件生成默认参数
struct RasterInfo *deafultRasterFromPoints(struct Points *points);
// 根据扇区类型搜索用于克里金插值的点
int searchNeighborhoods(double x, double y, int *neighbords, struct Points *points, struct NeighborhoodOption *neighborOpt, struct NeighbordSectorsWrap *sectorsWrap, struct DistanceAngleToCentroid *distanceAngleToCentroidCache);
// 区间分组：用于计算适配模型的参数
void binning(struct Points *points, struct SemivariogramOption *semivarOpt, int *lagsSemivarCount, double *lagsSemivarSum, struct PointPairDistance *pointPairDistance);
// 在qsort函数里边用于比较两个PointPixelDistance结构体数组角度的函数
int angleToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b);
// 在qsort函数里边用于比较两个PointPixelDistance结构体数组距离的函数
int distanceToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b);
// 根据临近点搜索条件生成sectors
struct NeighbordSectorsWrap *makeSectors(struct NeighborhoodOption *neighborOpt);

// struct PointPixelDistance **pointPixelDistance = NULL;

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

  // 2. 计算扇区参数
  struct DistanceAngleToCentroid *distanceAngleToCentroidCache = malloc(sizeof(struct DistanceAngleToCentroid) * pointNumbers);
  struct NeighbordSectorsWrap *sectorsWrap = makeSectors(neighborOpt);

  // 3. 栅格点插值
  int *neightbors = malloc(sizeof(int) * pointNumbers);
  double *rasterArray_ = malloc(sizeof(double) * rasterInfo->xsize * rasterInfo->ysize); // 栅格以行存储
  double *rasterArray = rasterArray_;

  for (int x = 0; x < rasterInfo->xsize; ++x)
  {
    double ptx = rasterInfo->left + rasterInfo->resolution * (x + 0.5);
    for (int y = 0; y < rasterInfo->ysize; ++y)
    {
      double pty = rasterInfo->top - rasterInfo->resolution * (y + 0.5);
      // 搜索用于插值的临近样本点
      int neighordsCount = searchNeighborhoods(ptx, pty, neightbors, points, neighborOpt, sectorsWrap, distanceAngleToCentroidCache);
      if (neighordsCount)
      {
        // TODO: 获取到临近点之后进行更进一步的计算
      }
    }
  }

  free(rasterArray_);
  free(distanceAngleToCentroidCache);
  free(sectorsWrap->sectors);
  free(sectorsWrap);
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

int searchNeighborhoods(double x, double y, int *neighbords, struct Points *points, struct NeighborhoodOption *neighborOpt, struct NeighbordSectorsWrap *sectorsWrap, struct DistanceAngleToCentroid *distanceAngleToCentroidCache)
{
  const double *pointsData = points->data;
  const int pointsNumber = points->numbers;
  const double maxDistance = neighborOpt->maxDistance < 0 ? DBL_MAX : neighborOpt->maxDistance;
  const int sectorType = neighborOpt->sectorType;
  const int minNeighbords = neighborOpt->minNeighbords;
  const int maxNeighbords = neighborOpt->maxNeighbords;

  // 首先计算栅格像素中心点（x, y）同各样本点的距离 & 角度
  int alternativeNeighbordsCount = 0;
  for (int i = 0; i < pointsNumber; ++i)
  {
    int pointBaseIndex = i * 3;
    // 样本点的x, y坐标
    double pointX = pointsData[pointBaseIndex], pointY = pointsData[pointBaseIndex + 1];
    double distance = sqrt(pow(x - pointX, 2) + pow(y - pointY, 2));
    if (distance < maxDistance)
    {
      struct DistanceAngleToCentroid *ppxd = &distanceAngleToCentroidCache[alternativeNeighbordsCount++];
      ppxd->index = i;
      ppxd->distance = distance;

      pointX -= x;
      pointY -= y;
      if (pointX == 0 && pointY == 0)
      {
        // TODO: 处理特殊情况，样本点与栅格像元中心点重合
        // 目前做法是忽略该点
        --alternativeNeighbordsCount;
        continue;
      }
      else
      {
        double angle = -atan2(pointY, pointX);
        angle += (angle < 0 ? M_2PI : 0);
        angle += M_PI_2;
        angle -= (angle > M_2PI ? M_2PI : 0);

        ppxd->angle = angle;
      }
    }
  }

  // 记录按照样本点与中心点角度排序
  struct DistanceAngleToCentroid *angleAscendingCache[alternativeNeighbordsCount];
  // 记录按照样本点与中心点距离排序
  struct DistanceAngleToCentroid *distanceAscendingCache[alternativeNeighbordsCount];

  for (int i = 0; i < alternativeNeighbordsCount; ++i)
  {
    angleAscendingCache[i] = &distanceAngleToCentroidCache[i];
  }
  qsort(angleAscendingCache, alternativeNeighbordsCount, sizeof(struct DistanceAngleToCentroid *), angleToCentroidCmp);

  // Note: 样本点与中心点的角度范围为[0, 2PI)
  double pointsMinAngle = ((struct DistanceAngleToCentroid *)(angleAscendingCache[0]))->angle,
         pointsMaxAngle = ((struct DistanceAngleToCentroid *)(angleAscendingCache[alternativeNeighbordsCount - 1]))->angle;

  // 迭代每一个扇区，然后将满足角度条件的点按距离由近及远排序。
  int neighbordsCount = 0;
  int lastAngleEndIndex = 0;
  for (int i = 0; i < sectorsWrap->count; ++i)
  {
    struct NeighbordSector *sector = &sectorsWrap->sectors[i];
    double sectorAngleFrom = sector->angleFrom;
    double sectorAngleTo = sector->angleTo;
    int sectorAlternativeNeighbordsCount = 0;

    // 迭代每一个样本点，如果该样本点的角度落入扇区，则将其记录
    // 之后将这些样本点按照距离排序，再按照距离远近和扇区能够容纳的样本点数量将其纳入扇区

    for (int j = 0; j < alternativeNeighbordsCount; ++j)
    {
      double angle = ((struct DistanceAngleToCentroid *)(angleAscendingCache[j]))->angle;
      double angleAdd2PI = angle + M_2PI;
      if ((angle >= sectorAngleFrom && angle < sectorAngleTo) || (angleAdd2PI >= sectorAngleFrom && angleAdd2PI < sectorAngleTo))
      {
        distanceAscendingCache[sectorAlternativeNeighbordsCount++] = angleAscendingCache[j];
      }
    }

    if (!sectorAlternativeNeighbordsCount)
    {
      continue;
    }

    qsort(distanceAscendingCache, sectorAlternativeNeighbordsCount, sizeof(struct DistanceAngleToCentroid *), distanceToCentroidCmp);

    for (int j = 0; j < sectorAlternativeNeighbordsCount; ++j)
    {
      struct DistanceAngleToCentroid *point = distanceAscendingCache[j];
      if (sector->count < sector->maxPointCount)
      {
        neighbords[neighbordsCount++] = point->index;
        ++sector->count;
      }
      else
      {
        break;
      }
    }

    // 检查扇区角度范围内是否有样本点
    // TODO: 已知样本点与中心点的角度范围为[0, 2PI)，扇区角度范围并不一定在[0, 2P])内，所以会存在问题
    if (sectorAngleFrom > pointsMaxAngle || sectorAngleTo <= pointsMinAngle)
    {
      // 该扇区角度范围内没有样本点
      continue;
    }
  }
  // TODO: 某个扇区的样本点数目不够的话将其丢弃？ ArcGIS的做法并没有将其丢弃

  return neighbordsCount;
}

int angleToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b)
{
  double angleA = (*a)->angle, angleB = (*b)->angle;

  if (angleA < angleB)
    return -1;
  if (angleA > angleB)
    return 1;
  return 0;
}

int distanceToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b)
{
  double distanceA = (*a)->distance, distanceB = (*b)->distance;
  if (distanceA < distanceB)
    return -1;
  if (distanceA > distanceB)
    return 1;
  return 0;
}

struct NeighbordSectorsWrap *makeSectors(struct NeighborhoodOption *neighborOpt)
{
  int sectorsCount = 0;
  double sectorInitOffsetAngle = 0;

  switch (neighborOpt->sectorType)
  {
  case SECTOR_TYPE_1:
    sectorsCount = 1;
    break;
  case SECTOR_TYPE_4S:
    sectorsCount = 4;
    break;
  case SECTOR_TYPE_4S_45D:
    sectorsCount = 4;
    sectorInitOffsetAngle = M_PI_4;
    break;
  case SECTOR_TYPE_8:
    sectorsCount = 8;
    break;
  }

  struct NeighbordSector *sectors = malloc(sizeof(struct NeighbordSector) * sectorsCount);

  const double ANGLE_PER_SECTOR = 2 * M_PI / sectorsCount;

  for (int i = 0; i < sectorsCount; ++i)
  {
    struct NeighbordSector *sector = &sectors[i];
    sector->id = i;
    sector->count = 0;
    sector->minPointCount = neighborOpt->minNeighbords;
    sector->maxPointCount = neighborOpt->maxNeighbords;

    if (i > 0)
    {
      struct NeighbordSector *lastSector = &sectors[i - 1];
      sector->angleFrom = lastSector->angleTo;
      sector->angleTo = lastSector->angleTo + ANGLE_PER_SECTOR;
    }
    else
    {
      sector->angleFrom = sectorInitOffsetAngle;
      sector->angleTo = sectorInitOffsetAngle + ANGLE_PER_SECTOR;
    }
  }

#ifdef MAKE_SECTORS_DEBUG

  puts("sectors:");
  for (int i = 0; i < sectorsCount; ++i)
  {
    struct NeighbordSector *sector = &sectors[i];
    printf("\tid=%d, count=%d, angleFrom=%f, angleTo=%f\n", sector->id, sector->count, sector->angleFrom, sector->angleTo);
  }

#endif

  struct NeighbordSectorsWrap *sectorsWrap = malloc(sizeof(struct NeighbordSectorsWrap));

  sectorsWrap->sectors = sectors;
  sectorsWrap->count = sectorsCount;

  return sectorsWrap;
}
