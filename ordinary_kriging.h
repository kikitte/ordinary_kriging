#include "rasterio.h"
#include "variogram.h"

#ifndef ORDINARY_KRIGING_H
#define ORDINARY_KRIGING_H

#define SECTOR_TYPE_1 "OneSector"
#define SECTOR_TYPE_4S "FourSector"
#define SECTOR_TYPE_4S_45D "FourSector45Offset"
#define SECTOR_TYPE_8 "EightSector"

struct Points
{
  double *data;
  int numbers;
};

// 存储点的index
struct PointsIndex
{
  int *data;
  int numbers;
};

// 记录扇区的id及其包含的临近点数量
struct Sector
{
  int minNeighbords;  // 扇区内容许的最大样本点数量(includsive)
  int maxNeighbords;  // 扇区内容许的最小样本点数量(includsive)
  double angleFrom;   // 起始角度(includsive: 如果样本点角度等于起始角度则将该点纳入扇区)
  double angleTo;     // 终止角度(excludsive: 如果样本点角度等于终止角度则不将该点纳入扇区)
  double maxDistance; // 扇区最大搜索距离
};

struct SectorsWrap
{
  struct Sector *sectors;
  int count;
};

// // 记录样本点到像元中心点的角度与距离
struct PointPairDistance
{
  int fromIndex;
  int toIndex;
  double distance;
};

// 记录像元中心点和插值样本点的距离与角度（连线与网格北方向的夹角）
struct DistanceAngle
{
  int index;
  double distance;
  double angle;
};

struct NeighborhoodOption
{
  char *SECTOR_TYPE;   // 临近点搜索的扇区类型: one_sector, four_sector, four_sector_45D_offset, eight_sector
  int MIN_NEIGHBORDS;  // 最小临近点数量
  int MAX_HEIGHBORDS;  // 最大临近点数量
  double MAX_DISTANCE; // 每个扇区最大搜索距离
};

void *ordinary_kriging(struct Points *points, struct VariogramModel *variogramOpt, struct SectorsWrap *sectorsWrap, struct RasterInfo *rasterInfo, VariogramFunction modelFunction);

#endif