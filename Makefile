OUTPUT_DIR = build

PROGRAM = ordinary_kriging

CC = gcc 

LDLIBS = -lm -lgdal -lpthread

# -g Produce debugging information in the operating system's native format (stabs, COFF, XCOFF, or DWARF).
CFLAGS = -std=c99 -O2

FILES = main.c config.c ordinary_kriging.c util.c variogram.c \
				math/gaussian.c math/mat_ops.c \
				cJSON/cJSON.c cJSON/cJSON_Utils.c

all:
	mkdir -p build
	${CC} ${CFLAGS} ${LDLIBS} ${FILES} -o ${OUTPUT_DIR}/${PROGRAM}
