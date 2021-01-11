#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"
#include "gaussian_jordan/gaussian.h"
#include "gaussian_jordan/mat_ops.h"

// #define MAKE_SECTORS_DEBUG

// 根据扇区类型搜索用于克里金插值的点
int searchNeighborhoods(double cellX, double cellY, int *neighbords, double *neighbordsDistance, struct Points *points, struct NeighborhoodOption *neighborOpt, struct NeighbordSectorsWrap *sectorsWrap, struct DistanceAngleToCentroid *distanceAngleToCentroidCache);
// 在qsort函数里边用于比较两个PointPixelDistance结构体数组角度的函数
int angleToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b);
// 在qsort函数里边用于比较两个PointPixelDistance结构体数组距离的函数
int distanceToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b);
// 根据临近点搜索条件生成sectors
struct NeighbordSectorsWrap *makeSectors(struct NeighborhoodOption *neighborOpt);

// struct PointPixelDistance **pointPixelDistance = NULL;

double *ordinaryKriging(struct Points *points,
                        struct VariogramModel *variogramOpt,
                        struct NeighborhoodOption *neighborOpt,
                        struct RasterInfo *rasterInfo)
{
  // 获取半变异函数和参数
  double (*modelFnc)(double H, double C0, double CX, double A) = NULL;
  switch (variogramOpt->VAR)
  {
  case VARIOGRAM_MODEL_LINEAR:
    return NULL;
  case VARIOGRAM_MODEL_SPHERICAL:
    modelFnc = sphericalVariogram;
    break;
  case VARIOGRAM_MODEL_EXPONENTIAL:
    modelFnc = exponentialVariogram;
    break;
  case VARIOGRAM_MODEL_GAUSSIAN:
    modelFnc = gaussianVariogram;
    break;
  default:
    return NULL;
  }

  // 2. 计算扇区参数
  struct DistanceAngleToCentroid *distanceAngleToCentroidCache = malloc(sizeof(struct DistanceAngleToCentroid) * points->numbers);
  struct NeighbordSectorsWrap *sectorsWrap = makeSectors(neighborOpt);

  // 3. 栅格点插值
  int *neightbors = malloc(sizeof(int) * points->numbers);
  double *neightborsDistance = malloc(sizeof(double) * points->numbers);
  double *rasterArray_ = malloc(sizeof(double) * rasterInfo->cols * rasterInfo->rows); // 栅格以行存储
  double *rasterArray = rasterArray_;

  for (int r = 0; r < rasterInfo->rows; ++r)
  {
    double cellY = rasterInfo->top - rasterInfo->resolution * (r + 0.5);
    for (int c = 0; c < rasterInfo->cols; ++c)
    {
      double cellX = rasterInfo->left + rasterInfo->resolution * (c + 0.5);
      // 搜索用于插值的临近样本点
      int neighordsCount = searchNeighborhoods(cellX, cellY, neightbors, neightborsDistance, points, neighborOpt, sectorsWrap, distanceAngleToCentroidCache);
      if (neighordsCount)
      {
        int dimen = neighordsCount + 1;

        int pointBaseIndex = 0;
        // AX = B, 求A的逆矩阵
        double *neighbordValues = malloc(sizeof(double) * neighordsCount);
        double **A = mat_zeros(dimen, dimen);
        double **B = mat_zeros(dimen, 1);

        for (int rr = 0; rr < dimen; ++rr)
        {
          A[rr][neighordsCount] = 1;
          A[rr][rr] = 0;
          B[rr][0] = rr < neighordsCount ? modelFnc(neightborsDistance[rr], variogramOpt->C0, variogramOpt->CX, variogramOpt->A) : 1;

          double nx, ny;
          if (rr < neighordsCount)
          {
            pointBaseIndex = neightbors[rr] * 3;
            nx = points->data[pointBaseIndex], ny = points->data[++pointBaseIndex];
            neighbordValues[rr] = points->data[++pointBaseIndex];
          }

          for (int cc = rr + 1; cc < neighordsCount; ++cc)
          {
            pointBaseIndex = neightbors[cc] * 3;
            double mx = points->data[pointBaseIndex], my = points->data[++pointBaseIndex];
            double mnx = mx - nx, mny = my - ny;
            double distance = sqrt(mnx * mnx + mny * mny);

            A[rr][cc] = modelFnc(distance, variogramOpt->C0, variogramOpt->CX, variogramOpt->A);
          }
        }

        // 将矩阵的右上部分复制到左下部分
        for (int rr = 0; rr < dimen; ++rr)
        {
          for (int cc = rr + 1; cc < dimen; ++cc)
          {
            A[cc][rr] = A[rr][cc];
          }
        }

        // 求解逆矩阵
        double **invA = inv(dimen, A);
        double **W = mat_mul(dimen, dimen, 1, invA, B);
        double cellValue = 0;
        for (int i = 0; i < neighordsCount; ++i)
        {
          cellValue += W[i][0] * neighbordValues[i];
        }
        *rasterArray++ = cellValue;

        printf("%f\n", cellValue);

        free(neighbordValues);
        free_ptr(dimen, A);
        free_ptr(dimen, B);
        free_ptr(dimen, W);
        free_ptr(dimen, invA);
      }
      else
      {
        *rasterArray++ = rasterInfo->nodata;
      }
      // 重置sectors计数
      for (int i = 0, iLen = sectorsWrap->count; i < iLen; ++i)
      {
        sectorsWrap->sectors[i].count = 0;
      }

    }
  }

  free(distanceAngleToCentroidCache);
  free(sectorsWrap->sectors);
  free(sectorsWrap);

  return rasterArray_;
}

int searchNeighborhoods(double cellX, double cellY, int *neighbords, double *neighbordsDistance, struct Points *points, struct NeighborhoodOption *neighborOpt, struct NeighbordSectorsWrap *sectorsWrap, struct DistanceAngleToCentroid *distanceAngleToCentroidCache)
{
  double *pointsData = points->data;
  const int pointsNumber = points->numbers;
  const double maxDistance = neighborOpt->maxDistance;
  const int sectorType = neighborOpt->sectorType;
  const int minNeighbords = neighborOpt->minNeighbords;
  const int maxNeighbords = neighborOpt->maxNeighbords;

  // 首先计算栅格像素中心点（x, y）同各样本点的距离 & 角度
  int alternativeNeighbordsCount = 0;
  for (int i = 0; i < pointsNumber; ++i, ++pointsData)
  {
    // 样本点相对于像元中心点的坐标
    double pointX = *pointsData++ - cellX, pointY = *pointsData++ - cellY;
    double distance = sqrt(pointX * pointX + pointY * pointY);
    // Note: ArcGIS好像超出范围了也要搜索如果扇区里点数不够的话
    // 这里不采用这种策略
    if (distance < maxDistance)
    {
      struct DistanceAngleToCentroid *ppxd = &distanceAngleToCentroidCache[alternativeNeighbordsCount++];
      ppxd->index = i;
      ppxd->distance = distance;

      if (distance)
      {
        double angle = -atan2(pointY, pointX);
        angle += (angle < 0 ? M_2PI : 0);
        angle += M_PI_2;
        angle -= (angle > M_2PI ? M_2PI : 0);

        ppxd->angle = angle;
      }
      else
      {
        // TODO: 处理特殊情况，样本点与栅格像元中心点重合
        // 目前做法是忽略该点
        --alternativeNeighbordsCount;
        continue;
      }
    }
  }

  // 记录按照样本点与中心点角度排序
  struct DistanceAngleToCentroid *angleAscendingCache[alternativeNeighbordsCount];
  // 记录按照样本点与中心点距离排序
  struct DistanceAngleToCentroid *distanceAscendingCache[alternativeNeighbordsCount];

  for (int i = 0; i < alternativeNeighbordsCount; ++i)
  {
    angleAscendingCache[i] = distanceAngleToCentroidCache + i;
  }
  qsort(angleAscendingCache, alternativeNeighbordsCount, sizeof(struct DistanceAngleToCentroid *), angleToCentroidCmp);

  // Note: 样本点与中心点的角度范围为[0, 2PI)
  double pointsMinAngle = ((struct DistanceAngleToCentroid *)(angleAscendingCache[0]))->angle,
         pointsMaxAngle = ((struct DistanceAngleToCentroid *)(angleAscendingCache[alternativeNeighbordsCount - 1]))->angle;

  // 迭代每一个扇区，然后将满足角度条件的点按距离由近及远排序。
  int neighbordsCount = 0;
  for (int i = 0, iLen = sectorsWrap->count; i < iLen; ++i)
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

    if (sectorAlternativeNeighbordsCount)
    {
      qsort(distanceAscendingCache, sectorAlternativeNeighbordsCount, sizeof(struct DistanceAngleToCentroid *), distanceToCentroidCmp);
      for (int j = 0; j < sectorAlternativeNeighbordsCount; ++j)
      {
        struct DistanceAngleToCentroid *point = distanceAscendingCache[j];
        if (sector->count < sector->maxPointCount)
        {
          neighbords[neighbordsCount] = point->index;
          neighbordsDistance[neighbordsCount] = point->distance;
          ++neighbordsCount;
          ++sector->count;
        }
        else
        {
          break;
        }
      }
    }
  }
  // TODO: 某个扇区的样本点数目不够的话将其丢弃？ ArcGIS的做法并没有将其丢弃

  return neighbordsCount;
}

int angleToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b)
{
  double angleA = (*a)->angle, angleB = (*b)->angle;

  return (angleA > angleB) - (angleA < angleB);
}

int distanceToCentroidCmp(struct DistanceAngleToCentroid **a, struct DistanceAngleToCentroid **b)
{
  double distanceA = (*a)->distance, distanceB = (*b)->distance;

  return (distanceA > distanceB) - (distanceA < distanceB);
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