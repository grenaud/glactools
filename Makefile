
CXX      = g++ -std=c++11 -fpermissive #-g -pg  #-g  # -fstack-protector 
#LIBGAB   = libgab/
#BAMTOOLS= $(realpath bamtools/)
#LIBGAB= $(realpath libgab/)


CXXFLAGS = -Wall -lm -O3 -lz -Ihtslib/ -Isamtools/ -Itabixpp/ -Ilibgab/ -Ilibgab/gzstream/ -Ibamtools/src/ -c #-Ibamtools/include/ -Ibamtools/src/
LDFLAGS  =   -lpthread -lm -lcurl -lbz2 -llzma -lz
LDLIBS   =     htslib/libhts.a samtools/libbam.a samtools/libst.a 

#LDFLAGS  =  ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o -lpthread -lm -lbz2 -llzma -lz


all: libgab/utils.o tabixpp/tabix.o samtools/bedidx.o bamtools/src/bamtools_fasta.o glactools 

%.o: %.cpp 
	${CXX} ${CXXFLAGS} $^ -o $@

libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

libgab/utils.o: htslib/libhts.so libgab/utils.h #bamtools/build/src/api/libbamtools.a 
	cd libgab/ && make utils.o && make -C gzstream/  && cd ..

tabixpp/tabix.hpp:
	rm -rf tabixpp/
	git clone --recursive https://github.com/grenaud/tabixpp.git

tabixpp/tabix.o: tabixpp/tabix.hpp htslib/libhts.so
	make -C tabixpp


bamtools/src/bamtools_fasta.o:
	cd bamtools/src/
	g++ -lm -c -I.  -Ishared/ -Iutils/ utils/bamtools_fasta.cpp
	cd ../..

#bamtools/src/api/BamAlignment.:
#	rm -rf bamtools/
#	git clone --recursive https://github.com/pezmaster31/bamtools.git

#bamtools/build/src/api/libbamtools.a: bamtools/src/api/BamAlignment.h
#	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..
#
#bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o: bamtools/build/src/api/libbamtools.a

samtools/bedidx.o:  samtools/sam.h
	cd samtools/ && make && cd ..

samtools/sam.h:
	rm -rf samtools/
	git clone --recursive https://github.com/samtools/samtools.git samtools/


htslib/libhts.so:  htslib/hts_internal.h
	cd htslib/ && make && cd ../

htslib/hts_internal.h:
	rm -rf htslib/
	git clone --recursive https://github.com/samtools/htslib.git

lib/libglactools.a: GlacIndex.o GlacIDXSTATS.o  RandomGenomicCoord.o GenomicWindows.o GlacWindows.o GlacWriter.o MultiVCFreader.o GLF2ACF.o VcfMulti2ACF.o T3andme2ACF.o AXT2ACF.o Vcf2ACF.o Vcf2GLF.o  SingleAllele.o BAM2ACF.o EIGENSTRAT2ACF.o BEAGLE2GLF.o BPLINK2ACF.o AlleleRecords.o SingleGL.o GlacUndef.o  GlacMindef.o GlacSegsite.o GlacRename.o GlacSharing.o GlacNoSharing.o GlacNoStrictSharing.o GlacNosingle.o GlacBedfilter.o Glac2FREQSPEC.o Glac2BED.o Glac2VCF.o ACF2BPLINK.o ACF2FASTA.o ACF2GPHOCS.o ACF2NEXUS.o ACF2TREEMIX.o ACF2GROSS.o ACF2EIGENSTRAT.o ACF2BETASCAN.o GlacViewer.o GlacParser.o GlacStats.o GlacReplaceAncestor.o GlacUsePopAsRootAnc.o GlacCAT.o VCFreader.o SimpleVCF.o CoreVCF.o ReadTabix.o SetVCFFilters.o GlacMeld.o GlacDown.o GlacPopsub.o GlacCompute.o SumStatD.o DstatResult.o DstatCounter.o Dstat_core.o SumStatF3.o F3Result.o F3Counter.o F3_core.o SumStatF2.o F2Result.o F2Counter.o F2_core.o SumStatAvgCoa.o AvgCoaAlleleCounter.o AvgCoaResult.o ComputeAvgCoa_core.o ComputeFst_core.o DistAlleleCounter.o DistResult.o SumStatDist.o ComputeDist_core.o FstAlleleCounter.o FstResult.o SumStatFst.o GlacClosest.o GlacRemovepop.o GlacIntersect.o GlacUnion.o GlacReheader.o FilterVCF.o GlactoolsOperations.o GenomicRange.o tabixpp/tabix.o libgab/utils.o bamtools/src/bamtools_fasta.o  samtools/libbam.a samtools/libst.a libgab/gzstream/gzstream.o #bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o
	ar rs lib/libglactools.a $^

glactools:	glactools.o  lib/libglactools.a htslib/libhts.a samtools/bedidx.o #bamtools/build/src/api/libbamtools.a
	${CXX} -o $@ $^ $(LDLIBS)  $(LDFLAGS)

clean :
	rm -f *.o glactools

.PHONY: test

test:	all
	cd test/ && bash testAll.sh && cd ..
