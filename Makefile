
CXX      = g++ -std=c++11
#LIBGAB   = libgab/
BAMTOOLS= $(realpath bamtools/)
LIBGAB= $(realpath libgab/)


CXXFLAGS = -Wall -lm -O3 -lz -I htslib/ -I tabixpp/ -I${LIBGAB} -I${LIBGAB}/gzstream/ -c
LDFLAGS  =  -lpthread -lm -lbz2 -llzma -lz


all: glactools 

%.o: %.cpp libgab/utils.o
	${CXX} ${CXXFLAGS} $^ -o $@

libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

libgab/utils.o: bamtools/lib/libbamtools.so htslib/libhts.so libgab/utils.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git

bamtools/lib/libbamtools.so: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

 htslib/libhts.so:  htslib/hts_internal.h
	cd htslib/ && make && cd ../

htslib/hts_internal.h:
	rm -rf htslib/
	git clone --recursive https://github.com/samtools/htslib.git


glactools:	glactools.o vcf2acf.o vcf2glf.o VCFreader.o SimpleVCF.o CoreVCF.o ReadTabix.o SetVCFFilters.o FilterVCF.o glactoolsOperations.o tabixpp/tabix.o htslib/libhts.a ${LIBGAB}/utils.o  libgab//gzstream/gzstream.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f *.o glactools

