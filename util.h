#include "ordinary_kriging.h"

#ifndef KRIGING_UTIL_h
#define KRIGING_UTIL_h

char *read_file_byte(const char *path);


CPLErr save_raster(struct RasterInfo *rastInfo, struct RasterInfosWrap *rasterInfosWrap, const char *output_filename);

#endif