
CXX      = g++ -std=c++11
#LIBGAB   = libgab/
BAMTOOLS= $(realpath bamtools/)
LIBGAB= $(realpath libgab/)


CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -c
LDFLAGS  = -lz


all: glactools 

libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

libgab/utils.o: bamtools/lib/libbamtools.so  libgab/utils.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git

bamtools/lib/libbamtools.so: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..


%.o: %.cpp ../lib/libgab/utils.o
	${CXX} ${CXXFLAGS} $^ -o $@

vcf2acf:	vcf2acf.o libgab/utils.o SetVCFFilters.o # tabix/libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o  libgab//gzstream/gzstream.o

glactools.o:	glactools.cpp
	${CXX} ${CXXFLAGS} glactools.cpp

glactools:	glactools.o vcf2acf.o SetVCFFilters.o ${LIBGAB}/utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f glactools.o glactools

