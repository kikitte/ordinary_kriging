#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ordinary_kriging.h"
#include "variogram.h"
#include "constants.h"
#include "math/gaussian.h"
#include "math/mat_ops.h"

/*
 * 根据扇区类型搜索用于克里金插值的点
 * 须知：对于扇区角度排序有依赖，扇区的角度范围必须随index顺序增长
 */
int search_neighbords(double cellX, double cellY, int *neighbords, double *neighbordsDistance, struct Points *points, struct SectorsWrap *sectorsWrap, struct DistanceAngle *distanceAngleCache, struct DistanceAngle **distanceAnglePointerCache);
int compare_angle(const void *, const void *);
void *convert_array(double *array, int arrayLength, GDALDataType gdt);

void *ordinary_kriging(void *param)
// void ordinary_kriging(void *);
{
  struct OrdinaryKrigingParams *ordinaryKrigingParam = param;
  struct Points *points = ordinaryKrigingParam->points;
  struct VariogramModel *variogramOpt = ordinaryKrigingParam->variogramOpt;
  struct SectorsWrap *sectorsWrap = ordinaryKrigingParam->sectorsWrap;
  struct RasterInfo *rasterInfo = ordinaryKrigingParam->rasterInfo;
  VariogramFunction modelFunction = ordinaryKrigingParam->modelFunction;

  // 1. 计算扇区参数
  struct DistanceAngle *distanceAngleCache = malloc(sizeof(struct DistanceAngle) * points->numbers);
  struct DistanceAngle **distanceAnglePointerCache = malloc(sizeof(struct DistanceAngle *) * points->numbers);

  // 2. 栅格点插值
  const int RASTER_COLS = rasterInfo->COLS, RASTER_ROWS = rasterInfo->ROWS;
  const int RASTER_CELL_NUMBERS = RASTER_COLS * RASTER_ROWS;

  int maxDimen = 2 + sectorsWrap->numbers * sectorsWrap->sectors[0].maxNeighbords;
  double **A = mat_zeros(maxDimen, maxDimen);    // 矩阵：系数矩阵和结果矩阵
  double *W = malloc(maxDimen * sizeof(double)); // 权重
  int *neighbords = malloc(sizeof(int) * points->numbers);
  double *neighbordsDistance = malloc(sizeof(double) * points->numbers);
  double *neighbordsValue = malloc(maxDimen * sizeof(double));
  double *rasterArray = malloc(sizeof(double) * RASTER_CELL_NUMBERS); // 栅格以行存储
  double *rasterArrayPtr = rasterArray;

  for (int r = 0; r < RASTER_ROWS; ++r)
  {
    double cellY = rasterInfo->TOPLEFT_Y - rasterInfo->RESOLUTION * (r + 0.5);
    for (int c = 0; c < RASTER_COLS; ++c)
    {
      double cellX = rasterInfo->TOPLEFT_X + rasterInfo->RESOLUTION * (c + 0.5);

      // 搜索用于插值的临近样本点
      int neighordsCount = search_neighbords(cellX, cellY, neighbords, neighbordsDistance, points, sectorsWrap, distanceAngleCache, distanceAnglePointerCache);
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
            A[rr][dimen] = modelFunction(neighbordsDistance[rr], variogramOpt->NUGGET, variogramOpt->SILL, variogramOpt->RANGE);

            // 记录点m和点m的坐标和x、y的差值和距离
            double nx, ny, mx, my, mnx, mny, distance;

            pointBaseIndex = neighbords[rr];
            nx = points->data[pointBaseIndex][0], ny = points->data[pointBaseIndex][1];
            neighbordsValue[rr] = points->data[pointBaseIndex][2];

            for (int cc = rr + 1; cc < neighordsCount; ++cc)
            {
              pointBaseIndex = neighbords[cc];
              mx = points->data[pointBaseIndex][0], my = points->data[pointBaseIndex][1];
              mnx = mx - nx, mny = my - ny;
              distance = sqrt(mnx * mnx + mny * mny);

              A[cc][rr] = A[rr][cc] = modelFunction(distance, variogramOpt->NUGGET, variogramOpt->SILL, variogramOpt->RANGE);
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
          *rasterArrayPtr++ = rasterInfo->NODATA_VALUE;
      }
      else
        *rasterArrayPtr++ = rasterInfo->NODATA_VALUE;
    }
  }

  free(W);
  free_ptr(maxDimen, A);
  free(distanceAngleCache);
  free(distanceAnglePointerCache);
  free(neighbords);
  free(neighbordsDistance);
  free(neighbordsValue);

  rasterArrayPtr = rasterArray;
  void *newRasterArray = convert_array(rasterArrayPtr, RASTER_CELL_NUMBERS, rasterInfo->GDAL_GDT);
  if (newRasterArray)
  {
    free(rasterArray);
    rasterInfo->rasterArray = newRasterArray;
  }
  else
  {
    rasterInfo->rasterArray = rasterArray;
  }

  return NULL;
}

int search_neighbords(double cellX, double cellY, int *neighbords, double *neighbordsDistance, struct Points *points, struct SectorsWrap *sectorsWrap, struct DistanceAngle *distanceAngleCache, struct DistanceAngle **distanceAnglePointerCache)
{
  const double(*pointsData)[3] = points->data;
  const int pointsNumber = points->numbers;

  // 首先计算栅格像素中心点（x, y）同各样本点的距离 & 角度
  for (int i = 0; i < pointsNumber; ++i)
  {
    // 样本点相对于像元中心点的坐标
    double pointX = pointsData[i][0] - cellX, pointY = pointsData[i][1] - cellY;
    double distance = sqrt(pointX * pointX + pointY * pointY);
    // 样本点与栅格像元中心点重合. atan2的两个参数均为0, 这时通常返回0， 所以这是可以接受的而无须特殊处理, 参考： https://stackoverflow.com/questions/47909048/what-will-be-atan2-output-for-both-x-and-y-as-0
    double angle = -atan2(pointY, pointX);
    // Note: 样本点与中心点的角度范围为[0, 2PI), 像元中心与样本点连线与网格北的夹角
    angle += ((angle >= -M_PI_2) ? M_PI_2 : (M_2PI + M_PI_2));

    struct DistanceAngle *pDistanceAngle = distanceAnglePointerCache[i] = distanceAngleCache + i;
    pDistanceAngle->index = i;
    pDistanceAngle->distance = distance;
    pDistanceAngle->angle = angle;
  }

  // 如果最后一个扇区的夹角大于2PI
  if (sectorsWrap->sectors[sectorsWrap->numbers - 1].angleTo > M_2PI)
  {
    const double firstSectorAngleFrom = sectorsWrap->sectors[0].angleFrom;
    for (int i = 0; i < pointsNumber; ++i)
    {
      if (distanceAnglePointerCache[i]->angle < firstSectorAngleFrom)
      {
        distanceAnglePointerCache[i]->angle += M_2PI;
      }
    }
  }
  // 将样本点按照夹角排序
  qsort(distanceAnglePointerCache, pointsNumber, sizeof(struct DistanceAngle *), compare_angle);

  // 迭代每一个扇区，然后将满足角度条件的点按距离由近及远排序。
  int neighbordsCount = 0, firstPointIndex = 0, lastPointIndex = 0;
  for (int i = 0; i < sectorsWrap->numbers; ++i)
  {
    const struct Sector *sector = &sectorsWrap->sectors[i];
    const double sectorAngleFrom = sector->angleFrom;
    const double sectorAngleTo = sector->angleTo;
    const int sectorMinNeighbors = sector->minNeighbords;
    const int sectorMaxNeighbords = sector->maxNeighbords;
    int sectorNeighbords = 0;

    firstPointIndex = lastPointIndex;
    while (lastPointIndex < pointsNumber)
    {
      double angle = distanceAnglePointerCache[lastPointIndex]->angle;
      if (angle >= sectorAngleFrom && angle < sectorAngleTo)
      {
        ++lastPointIndex;
      }
      else
      {
        break;
      }
    }

    // 满足条件的样本点的计数
    int sectorPointNumbers = lastPointIndex - firstPointIndex;
    int maxSelectionNumbers = sectorMaxNeighbords < sectorPointNumbers ? sectorMaxNeighbords : sectorPointNumbers;
    int minSelectionNumbers = sectorMinNeighbors < sectorPointNumbers ? sectorMinNeighbors : sectorPointNumbers;

    // 记录按照样本点与中心点距离排序
    for (int j = firstPointIndex, maxJ = j + maxSelectionNumbers; j < maxJ; ++j)
    {
      // 使用选择排序找到距离最近的样本点
      int minDistanceAnglePtrIndex = j;
      struct DistanceAngle *minDistanceAnglePtr = distanceAnglePointerCache[j];
      for (int k = j + 1; k < lastPointIndex; ++k)
      {
        if (distanceAnglePointerCache[k]->distance < minDistanceAnglePtr->distance)
        {
          minDistanceAnglePtrIndex = k;
          minDistanceAnglePtr = distanceAnglePointerCache[k];
        }
      }
      distanceAnglePointerCache[minDistanceAnglePtrIndex] = distanceAnglePointerCache[j];
      distanceAnglePointerCache[j] = minDistanceAnglePtr;

      if (minDistanceAnglePtr->distance > sector->maxDistance)
      {
        if (sectorNeighbords >= sectorMinNeighbors)
        {
          break;
        }
        else
        {
          maxJ = firstPointIndex + minSelectionNumbers;
        }
      }

      neighbords[neighbordsCount] = minDistanceAnglePtr->index;
      neighbordsDistance[neighbordsCount] = minDistanceAnglePtr->distance;
      ++neighbordsCount;
      ++sectorNeighbords;
    }
  }
  // TODO: 某个扇区的样本点数目不够的话将其丢弃？ ArcGIS的做法并没有将其丢弃

  return neighbordsCount;
}

int compare_angle(const void *a, const void *b)
{
  double angleA = (*(struct DistanceAngle **)a)->angle;
  double angleB = (*(struct DistanceAngle **)b)->angle;

  return (angleA > angleB) - (angleA < angleB);
}

void *convert_array(double *array, int arrayLength, GDALDataType gdt)
{
  switch (gdt)
  {
  case GDT_Byte:
  {
    char *rasterArrayChar = malloc(sizeof(char) * arrayLength);
    char *rasterArrayCharPtr = rasterArrayChar;
    for (int i = 0; i < arrayLength; ++i)
      *rasterArrayCharPtr++ = *array++;
    return rasterArrayChar;
  }
  case GDT_Int16:
  {
    int16_t *rasterArrayInt16 = malloc(sizeof(int16_t) * arrayLength);
    int16_t *rasterArrayInt16Ptr = rasterArrayInt16;
    for (int i = 0; i < arrayLength; ++i)
      *rasterArrayInt16Ptr++ = *array++;
    return rasterArrayInt16;
  }
  case GDT_Int32:
  {
    int32_t *rasterArrayInt32 = malloc(sizeof(int32_t) * arrayLength);
    int32_t *rasterArrayInt32Ptr = rasterArrayInt32;
    for (int i = 0; i < arrayLength; ++i)
      *rasterArrayInt32Ptr++ = *array++;
    return rasterArrayInt32;
  }
  case GDT_UInt16:
  {
    uint16_t *rasterArrayUInt16 = malloc(sizeof(uint16_t) * arrayLength);
    uint16_t *rasterArrayUInt16Ptr = rasterArrayUInt16;
    for (int i = 0; i < arrayLength; ++i)
      *rasterArrayUInt16Ptr++ = *array++;
    return rasterArrayUInt16;
  }
  case GDT_UInt32:
  {
    uint32_t *rasterArrayUInt32 = malloc(sizeof(uint32_t) * arrayLength);
    uint32_t *rasterArrayUInt32Ptr = rasterArrayUInt32;
    for (int i = 0; i < arrayLength; ++i)
      *rasterArrayUInt32Ptr++ = *array++;
    return rasterArrayUInt32;
  }
  case GDT_Float32:
  {
    float *rasterArrayFloat = malloc(sizeof(float) * arrayLength);
    float *rasterArrayFloatPtr = rasterArrayFloat;
    for (int i = 0; i < arrayLength; ++i)
      *rasterArrayFloatPtr++ = *array++;
    return rasterArrayFloat;
  }
  default:
    return NULL;
  }
}