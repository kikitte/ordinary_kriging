#include "cJSON/cJSON.h"
#include "ordinary_kriging.h"
#include "variogram.h"

#ifndef KRIGING_CONFIG
#define KRIGING_CONFIG

// 解析配置文件
cJSON *parse_config(const char *configPath);
// 解析VariogramModel
struct VariogramModel *get_variogram(cJSON *config);
// 解析NeighborhoodOption
struct NeighborhoodOption *get_neighborhood(cJSON *config);
// 解析RasterInfo
struct RasterInfo *get_raster(cJSON *config);
// 获取输入样本点的路径
char *get_input_points_path(cJSON *config);
// 或缺输出栅格的路径
char *get_output_raster_path(cJSON *config);
#endif