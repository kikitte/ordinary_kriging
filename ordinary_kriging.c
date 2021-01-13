#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"
#include "math/gaussian.h"
#include "math/mat_ops.h"

// #define MAKE_SECTORS_DEBUG

// 根据扇区类型搜索用于克里金插值的点
int searchNeighborhoods(double cellX, double cellY, int *neighbords, double *neighbordsDistance, struct Points *points, struct NeighborhoodOption *neighborOpt, struct SectorsWrap *sectorsWrap, struct DistanceAngle *distanceAngleCache, struct DistanceAngle **distanceAnglePointerCache);
// 在qsort函数里边用于比较两个PointPixelDistance结构体数组距离的函数
int distanceToCentroidCmp(const void *a, const void *b);
// 根据临近点搜索条件生成sectors
struct SectorsWrap *makeSectors(struct NeighborhoodOption *neighborOpt);

double *ordinaryKriging(struct Points *points,
                        struct VariogramModel *variogramOpt,
                        struct NeighborhoodOption *neighborOpt,
                        struct RasterInfo *rasterInfo)
{
  // 获取半变异函数
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
  struct SectorsWrap *sectorsWrap = makeSectors(neighborOpt);
  struct DistanceAngle *distanceAngleCache = malloc(sizeof(struct DistanceAngle) * points->numbers);
  struct DistanceAngle **distanceAnglePointerCache = malloc(sizeof(struct DistanceAngle *) * points->numbers);

  // 3. 栅格点插值
  int maxDimen = 2 + sectorsWrap->count * sectorsWrap->sectors[0].maxPointCount;
  double **A = mat_zeros(maxDimen, maxDimen);    // 矩阵：系数矩阵和结果矩阵
  double *W = malloc(maxDimen * sizeof(double)); // 权重
  int *neighbords = malloc(sizeof(int) * points->numbers);
  double *neighbordsDistance = malloc(sizeof(double) * points->numbers);
  double *neighbordsValue = malloc(maxDimen * sizeof(double));
  double *rasterArray = malloc(sizeof(double) * rasterInfo->cols * rasterInfo->rows); // 栅格以行存储
  double *rasterArrayPtr = rasterArray;

  for (int r = 0; r < rasterInfo->rows; ++r)
  {
    double cellY = rasterInfo->top - rasterInfo->resolution * (r + 0.5);
    for (int c = 0; c < rasterInfo->cols; ++c)
    {
      double cellX = rasterInfo->left + rasterInfo->resolution * (c + 0.5);

      // 搜索用于插值的临近样本点
      int neighordsCount = searchNeighborhoods(cellX, cellY, neighbords, neighbordsDistance, points, neighborOpt, sectorsWrap, distanceAngleCache, distanceAnglePointerCache);
      if (neighordsCount)
      {
        int pointBaseIndex = 0;
        int dimen = neighordsCount + 1;

        // 生成用于计算权重的矩阵
        for (int rr = 0; rr < dimen; ++rr)
        {
          A[neighordsCount][rr] = A[rr][neighordsCount] = 1;
          A[rr][rr] = 0;

          if (rr < neighordsCount)
          {
            A[rr][dimen] = modelFnc(neighbordsDistance[rr], variogramOpt->C0, variogramOpt->CX, variogramOpt->A);

            // 记录点m和点m的坐标和x、y的差值和距离
            double nx, ny, mx, my, mnx, mny, distance;

            pointBaseIndex = neighbords[rr] * 3;
            nx = points->data[pointBaseIndex], ny = points->data[++pointBaseIndex];
            neighbordsValue[rr] = points->data[++pointBaseIndex];

            for (int cc = rr + 1; cc < neighordsCount; ++cc)
            {
              pointBaseIndex = neighbords[cc] * 3;
              mx = points->data[pointBaseIndex], my = points->data[++pointBaseIndex];
              mnx = mx - nx, mny = my - ny;
              distance = sqrt(mnx * mnx + mny * mny);

              A[cc][rr] = A[rr][cc] = modelFnc(distance, variogramOpt->C0, variogramOpt->CX, variogramOpt->A);
            }
          }
          else
          {
            A[neighordsCount][dimen] = 1;
          }
        }

        // 求解
        if (gaussian_elimination_solve(dimen, A, W, 0.0001))
        {
          // 有解
          double cellValue = 0;
          for (int i = 0; i < neighordsCount; ++i)
            cellValue += W[i] * neighbordsValue[i];

          *rasterArrayPtr++ = cellValue;
        }
        else
          // 无解
          *rasterArrayPtr++ = rasterInfo->nodata;
      }
      else
        *rasterArrayPtr++ = rasterInfo->nodata;
    }
  }

  free(W);
  free_ptr(maxDimen, A);
  free(distanceAngleCache);
  free(distanceAnglePointerCache);
  free(sectorsWrap->sectors);
  free(sectorsWrap);
  free(neighbords);
  free(neighbordsDistance);
  free(neighbordsValue);

  return rasterArray;
}

int searchNeighborhoods(double cellX, double cellY, int *neighbords, double *neighbordsDistance, struct Points *points, struct NeighborhoodOption *neighborOpt, struct SectorsWrap *sectorsWrap, struct DistanceAngle *distanceAngleCache, struct DistanceAngle **distanceAnglePointerCache)
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
    // Note: ArcGIS好像超出范围了也要搜索如果扇区里点数不够的话, 这里不采用这种策略
    if (distance < maxDistance)
    {
      struct DistanceAngle *ppxd = distanceAngleCache + alternativeNeighbordsCount++;
      ppxd->index = i;
      ppxd->distance = distance;
      // 样本点与栅格像元中心点重合. atan2的两个参数均为0, 这时通常返回0， 所以这是可以接受的而无须特殊处理, 参考： https://stackoverflow.com/questions/47909048/what-will-be-atan2-output-for-both-x-and-y-as-0
      double angle = -atan2(pointY, pointX);
      // Note: 样本点与中心点的角度范围为[0, 2PI), 像元中心与样本点连线与网格北的夹角
      angle += ((angle >= -M_PI_2) ? M_PI_2 : (M_2PI + M_PI_2));

      ppxd->angle = angle;
    }
  }

  // 迭代每一个扇区，然后将满足角度条件的点按距离由近及远排序。
  int neighbordsCount = 0;
  for (int i = 0; i < sectorsWrap->count; ++i)
  {
    struct Sector *sector = &sectorsWrap->sectors[i];
    double sectorAngleFrom = sector->angleFrom;
    double sectorAngleTo = sector->angleTo;
    int sectorAlternativeNeighbordsCount = 0;

    // 迭代每一个样本点，如果该样本点的角度落入扇区，则将其记录
    // 之后将这些样本点按照距离排序，再按照距离远近和扇区能够容纳的样本点数量将其纳入扇区

    for (int j = 0; j < alternativeNeighbordsCount; ++j)
    {
      double angle = (distanceAngleCache + j)->angle;
      double angleAdd2PI = angle + M_2PI;
      if ((angle >= sectorAngleFrom && angle < sectorAngleTo) || (angleAdd2PI >= sectorAngleFrom && angleAdd2PI < sectorAngleTo))
        distanceAnglePointerCache[sectorAlternativeNeighbordsCount++] = distanceAngleCache + j;
    }

    if (sectorAlternativeNeighbordsCount)
    {
      // 满足条件的样本点的计数
      int count = 0, maxCount = sector->maxPointCount;
      // 记录按照样本点与中心点距离排序
      qsort(distanceAnglePointerCache, sectorAlternativeNeighbordsCount, sizeof(struct DistanceAngle *), distanceToCentroidCmp);
      for (int j = 0; j < sectorAlternativeNeighbordsCount; ++j)
      {
        struct DistanceAngle *point = distanceAnglePointerCache[j];
        if (count < maxCount)
        {
          neighbords[neighbordsCount] = point->index;
          neighbordsDistance[neighbordsCount] = point->distance;
          ++neighbordsCount;
          ++count;
        }
        else
          break;
      }
    }
  }
  // TODO: 某个扇区的样本点数目不够的话将其丢弃？ ArcGIS的做法并没有将其丢弃

  return neighbordsCount;
}

int distanceToCentroidCmp(const void *a, const void *b)
{
  double distanceA = (*(struct DistanceAngle **)a)->distance, distanceB = (*(struct DistanceAngle **)b)->distance;

  return (distanceA > distanceB) - (distanceA < distanceB);
}

struct SectorsWrap *makeSectors(struct NeighborhoodOption *neighborOpt)
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
  default:
    return NULL;
  }

  struct Sector *sectors = malloc(sizeof(struct Sector) * sectorsCount);

  const double ANGLE_PER_SECTOR = 2 * M_PI / sectorsCount;

  for (int i = 0; i < sectorsCount; ++i)
  {
    struct Sector *sector = &sectors[i];
    sector->id = i;
    sector->minPointCount = neighborOpt->minNeighbords;
    sector->maxPointCount = neighborOpt->maxNeighbords;

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

#ifdef MAKE_SECTORS_DEBUG

  puts("sectors:");
  for (int i = 0; i < sectorsCount; ++i)
  {
    struct NeighbordSector *sector = &sectors[i];
    printf("\tid=%d, count=%d, angleFrom=%f, angleTo=%f\n", sector->id, sector->count, sector->angleFrom, sector->angleTo);
  }

#endif

  struct SectorsWrap *sectorsWrap = malloc(sizeof(struct SectorsWrap));

  sectorsWrap->sectors = sectors;
  sectorsWrap->count = sectorsCount;

  return sectorsWrap;
}