
CXX      = g++ -std=c++11  -g -pg 
#LIBGAB   = libgab/
BAMTOOLS= $(realpath bamtools/)
LIBGAB= $(realpath libgab/)


CXXFLAGS = -Wall -lm -O3 -lz -Ihtslib/ -Ibamtools/include/ -Ibamtools/src/ -Itabixpp/ -I${LIBGAB} -I${LIBGAB}/gzstream/ -c
LDFLAGS  =  -lpthread -lm -lbz2 -llzma -lz


all: libgab/utils.o tabixpp/tabix.o glactools 

%.o: %.cpp 
	${CXX} ${CXXFLAGS} $^ -o $@

libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

libgab/utils.o: bamtools/lib/libbamtools.so htslib/libhts.so libgab/utils.h
	make -C libgab

tabixpp/tabix.hpp:
	rm -rf tabixpp/
	git clone --recursive https://github.com/grenaud/tabixpp.git

tabixpp/tabix.o: tabixpp/tabix.hpp htslib/libhts.so
	make -C tabixpp


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


glactools:	glactools.o GlacIndex.o GlacIDXSTATS.o GlacWriter.o MultiVCFreader.o GLF2ACF.o VcfMulti2ACF.o T3andme2ACF.o Vcf2ACF.o Vcf2GLF.o AlleleRecords.o SingleAllele.o BAM2ACF.o SingleGL.o GlacUndef.o GlacSegsite.o GlacSharing.o GlacNoSharing.o GlacNoStrictSharing.o GlacBedfilter.o Glac2FREQSPEC.o Glac2BED.o DstatResult.o DstatCounter.o Dstat_core.o ACF2BPLINK.o ACF2FASTA.o ACF2GPHOCS.o ACF2NEXUS.o ACF2TREEMIX.o ACF2EIGENSTRAT.o GlacViewer.o GlacParser.o GlacStats.o GlacReplaceAncestor.o GlacUsePopAsRootAnc.o GlacCAT.o VCFreader.o SimpleVCF.o CoreVCF.o ReadTabix.o SetVCFFilters.o GlacMeld.o GlacPopsub.o GlacCompute.o SumStatD.o SumStatAvgCoa.o AlleleCounter.o AvgCoaResult.o ComputeAvgCoa_core.o GlacClosest.o GlacRemovepop.o GlacIntersect.o GlacUnion.o GlacReheader.o FilterVCF.o GlactoolsOperations.o GenomicRange.o tabixpp/tabix.o htslib/libhts.a ${LIBGAB}/utils.o bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o bamtools/lib/libbamtools.a libgab/gzstream/gzstream.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f *.o glactools
