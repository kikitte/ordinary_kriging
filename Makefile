OUTPUT_DIR = build

PROGRAM = ordinary_kriging

CC = gcc 

LDLIBS = -lm -lgdal

# -g Produce debugging information in the operating system's native format (stabs, COFF, XCOFF, or DWARF).
CFLAGS = -std=c99 -O2

FILES = main.c config.c ordinary_kriging.c rasterio.c variogram.c \
				math/gaussian.c math/mat_ops.c \
				cJSON/cJSON.c cJSON/cJSON_Utils.c

all:
	${CC} ${CFLAGS} ${LDLIBS} ${FILES} -o ${OUTPUT_DIR}/${PROGRAM}
