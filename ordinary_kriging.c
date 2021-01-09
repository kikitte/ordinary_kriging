#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"

#define BINNING_DEBUG

// 在没有制定输出栅格参数时根据已有条件生成默认参数
struct RasterInfo *deafultRasterFromPoints(struct Points *points);
// 根据扇区类型搜索用于克里金插值的点
int searchNeighborhoods(double x, double y, int *neighbords, struct Points *points, struct NeighborhoodOption *neighborOpt, struct NeighbordSectorsWrap *sectorsWrap, struct PointPixelDistance *pointPixelDistanceCache);
// 区间分组：用于计算适配模型的参数
void binning(struct Points *points, struct SemivariogramOption *semivarOpt, int *lagsSemivarCount, double *lagsSemivarSum, struct PointPairDistance *pointPairDistance);
// 在qsort函数里边用于比较两个PointPixelDistance结构体数组距离的函数
int pointPixelDistanceCmp(struct PointPixelDistance *a, struct PointPixelDistance *b);
// 根据临近点搜索条件生成sectors
struct NeighbordSectorsWrap *makeSectors(struct NeighborhoodOption *neighborOpt);
// 将样本点分配到扇区，成功返回真值、失败返回假值
int distributePointIntoSectors(double distance, double angle, struct NeighbordSectorsWrap *sectorsWrap);

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
  // 即使pointPixelDistance被修改了，由于保存最原始的指针地址，我们也还是能正确释放他
  struct PointPixelDistance *pointPixelDistance = malloc(sizeof(struct PointPixelDistance) * pointNumbers);
  struct NeighbordSectorsWrap *sectorsWrap = makeSectors(neighborOpt);

  // 3. 栅格点插值
  int *neightborsIndex = malloc(sizeof(int) * pointNumbers);
  double *rasterArray_ = malloc(sizeof(double) * rasterInfo->xsize * rasterInfo->ysize); // 栅格以行存储
  double *rasterArray = rasterArray_;
  for (int x = 0; x < rasterInfo->xsize; ++x)
  {
    double ptx = rasterInfo->left + rasterInfo->resolution * (x + 0.5);
    for (int y = 0; y < rasterInfo->ysize; ++y)
    {
      double pty = rasterInfo->top - rasterInfo->resolution * (y + 0.5);
      // 搜索用于插值的临近样本点
      int neighordsCount = searchNeighborhoods(ptx, pty, neightborsIndex, points, neighborOpt, sectorsWrap, pointPixelDistance);
      if (neighordsCount)
      {
        printf("%d, %f, %f; ", neighordsCount, ptx, pty);
        // *rasterArray++ = result;
      }
    }
  }

  free(rasterArray_);
  free(pointPixelDistance);
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

int searchNeighborhoods(double x, double y, int *neighbords, struct Points *points, struct NeighborhoodOption *neighborOpt, struct NeighbordSectorsWrap *sectorsWrap, struct PointPixelDistance *pointPixelDistanceCache)
{
  // 重置pointPixelDistanceCache数组
  int alternativeNeighbordsCount;
  for (int i = 0; i < points->numbers; ++i)
  {
    pointPixelDistanceCache[i].index = -1;
  }

  const double *pointsData = points->data;
  const double maxDistance = neighborOpt->maxDistance < 0 ? DBL_MAX : neighborOpt->maxDistance;
  const int sectorType = neighborOpt->sectorType;
  const int minNeighbords = neighborOpt->minNeighbords;
  const int maxNeighbords = neighborOpt->maxNeighbords;

  // 首先计算栅格像素中心点（x, y）同各样本点的距离
  alternativeNeighbordsCount = 0;
  for (int i = 0; i < points->numbers; ++i)
  {
    int pointBaseIndex = i * 3;
    // 样本点的x, y坐标
    double pointX = pointsData[pointBaseIndex], pointY = pointsData[pointBaseIndex + 1];
    double distance = sqrt(pow(x - pointX, 2) + pow(y - pointY, 2));
    if (distance < maxDistance)
    {
      struct PointPixelDistance *ppxd = &pointPixelDistanceCache[alternativeNeighbordsCount++];
      ppxd->index = i;
      ppxd->distance = distance;
    }
  }
  // 到这里已经记录下所有满足距离条件的样本点index, 下一步根据距离进行对pointPixelDistance数组进行排序，得到距离从小到大的临近点排行
  struct PointPixelDistance **distanceAscendingCache[points->numbers];
  for (int i = 0; i < points->numbers; ++i)
  {
    distanceAscendingCache[i] = &pointPixelDistanceCache[i];
  }
  qsort(distanceAscendingCache, alternativeNeighbordsCount, sizeof(struct PointPixelDistance *), pointPixelDistanceCmp);
  // 从距离由近及远迭代alternativePoints，对每一个临近点尝试将其分入扇区，直到碰到不能分入任何扇区的点才结束
  int neighbordsCount = 0;
  for (int i = 0; i < alternativeNeighbordsCount; ++i)
  {
    struct PointPixelDistance *pointPixelDistance = distanceAscendingCache[i];
    int pointBaseIndex = pointPixelDistance->index * 3;
    double pointX = pointsData[pointBaseIndex], pointY = pointsData[pointBaseIndex + 1];
    double distance = pointPixelDistance->distance;
    double angle = asin(pointX / distance); // 样本点与像元中心点连线与网格北方向的夹角
    // 以下对角度进行一些转换以确保是与网格北方向的夹角
    if (angle > 0)
    {
      if (pointY < 0)
        angle += M_PI_2;
    }
    else if (angle < 0)
    {
      angle += M_2PI;
      if (pointY < 0)
        angle -= M_PI_2;
    }
    else
    {
      if (pointY > 0)
        angle = 0;
      else if (pointY < 0)
        angle = M_PI;
      else
      {
        // TODO： 特殊情况已知点与待插值点空间位置完全一样
        continue;
      }
    }

    if (distributePointIntoSectors(distance, angle, sectorsWrap))
    {
      neighbords[neighbordsCount++] = pointPixelDistance->index;
    }
    else
    {
      break;
    }
  }
  // TODO: 某个扇区的样本点数目不够的话将其丢弃？

  return neighbordsCount;
}

int distributePointIntoSectors(double distance, double angle, struct NeighbordSectorsWrap *sectorsWrap)
{
  for (int i = 0; i < sectorsWrap->count; ++i)
  {
    struct NeighbordSector *sector = &sectorsWrap->sectors[i];

    if ((angle >= sector->angleFrom && angle < sector->angleTo) || (angle + M_2PI >= sector->angleFrom && angle + M_2PI < sector->angleTo))
    {

      // 已经满足角度条件
      // 检查该扇区内样本点是否够数了
      if (sector->count < sector->maxPointCount)
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
  }
  return 0;
}

int pointPixelDistanceCmp(struct PointPixelDistance *a, struct PointPixelDistance *b)
{
  double aD = a->distance, bD = b->distance;
  return (aD > bD) - (aD < bD);
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

  for (int i = 0; i < sectorsCount; ++i)
  {
    struct NeighbordSector *sector = &sectors[i];
    printf("id=%d, count=%d, angleFrom=%f, angleTo=%f\n", sector->id, sector->count, sector->angleFrom, sector->angleTo);
  }

  struct NeighbordSectorsWrap *sectorsWrap = malloc(sizeof(struct NeighbordSectorsWrap));

  sectorsWrap->sectors = sectors;
  sectorsWrap->count = sectorsCount;

  return sectorsWrap;
}
