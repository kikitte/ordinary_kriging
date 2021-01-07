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

// 记录两点的index及其距离
struct PointPairDistance
{
  int fromIndex;
  int toIndex;
  double distance;
};

struct RasterInfo
{
  double top;
  double left;
  int xsize;
  int ysize;
  double resolution;
};

struct SemivariogramOption
{
  double lag;
  double lagNumbers;
  int modelType;
};

struct NeighborhoodOption
{
  int sectorType;     // 临近点搜索的扇区类型
  int minNeighbords;  // 最小临近点数量
  int maxNeighbords;  // 最大临近点数量
  double maxDistance; // 每个扇区最大搜索距离，负数和零表示无限大
};

double *ordinaryKriging(struct Points *points,
                        struct SemivariogramOption *semivarOpt,
                        struct NeighborhoodOption *neighborOpt,
                        struct RasterInfo *rasterInfo);