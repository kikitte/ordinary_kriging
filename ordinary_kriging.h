#include "variogram.h"

#define SECTOR_TYPE_1 5
#define SECTOR_TYPE_4S 6
// arcgis: 4 sectors with 45 degrees offset
#define SECTOR_TYPE_4S_45D 7
#define SECTOR_TYPE_8 8

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
struct NeighbordSector
{
  int id;
  int count;
  int minPointCount; // 扇区内容许的最大样本点数量(includsive)
  int maxPointCount; // 扇区内容许的最小样本点数量(includsive)
  double angleFrom;  // 起始角度(includsive: 如果样本点角度等于起始角度则将该点纳入扇区)
  double angleTo;    // 终止角度(excludsive: 如果样本点角度等于终止角度则不将该点纳入扇区)
};

struct NeighbordSectorsWrap
{
  struct NeighbordSector *sectors;
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
struct DistanceAngleToCentroid
{
  int index;
  double distance;
  double angle;
};

struct RasterInfo
{
  double top;
  double left;
  int cols;
  int rows;
  double resolution;
  const double nodata;
};

struct NeighborhoodOption
{
  int sectorType;     // 临近点搜索的扇区类型
  int minNeighbords;  // 最小临近点数量
  int maxNeighbords;  // 最大临近点数量
  double maxDistance; // 每个扇区最大搜索距离
};

double *ordinaryKriging(struct Points *points,
                        struct VariogramModel *semivarOpt,
                        struct NeighborhoodOption *neighborOpt,
                        struct RasterInfo *rasterInfo);