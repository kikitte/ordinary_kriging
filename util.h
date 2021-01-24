#include "ordinary_kriging.h"

#ifndef KRIGING_UTIL_h
#define KRIGING_UTIL_h

char *read_file_byte(const char *path);


CPLErr save_raster(void *raster_array, struct RasterInfo *rast_param, const char *output_filename);

#endif