
CXX      = g++   
LIBGLACTOOLS   = ../../../

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGLACTOOLS}/htslib/ -I${LIBGLACTOOLS}/bamtools/include/ -I${LIBGLACTOOLS}/bamtools/src/ -I${LIBGLACTOOLS}/tabixpp/ -I${LIBGLACTOOLS}/ -I${LIBGLACTOOLS}/libgab/ -I${LIBGLACTOOLS}/libgab/gzstream/  -c
LDFLAGS  = -lpthread -lm -lbz2 -llzma -lz


all: readACFGLFFile 

readACFGLFFile.o:	readACFGLFFile.cpp
	${CXX} ${CXXFLAGS} readACFGLFFile.cpp


readACFGLFFile:	readACFGLFFile.o ${LIBGLACTOOLS}/lib/libglactools.a ${LIBGLACTOOLS}/htslib/libhts.a ${LIBGLACTOOLS}/bamtools/lib/libbamtools.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f readACFGLFFile.o readACFGLFFile

