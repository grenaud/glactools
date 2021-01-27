
CXX      = g++ -std=c++11 -fpermissive #-g -pg  #-g  # -fstack-protector 
#LIBGAB   = libgab/
#BAMTOOLS= $(realpath bamtools/)
#LIBGAB= $(realpath libgab/)


CXXFLAGS = -Wall -Wno-unused-variable -Wno-unused-but-set-variable -lm -O3 -lz -Ihtslib/ -Isamtools/ -Itabixpp/ -Ilibgab/ -Ilibgab/gzstream/ -Ibamtools/src/ -c #-Ibamtools/include/ -Ibamtools/src/
LDFLAGS  =   -lpthread -lm -lcurl -lbz2 -llzma -lz
LDLIBS   =     htslib/libhts.a samtools/libbam.a samtools/libst.a 

#LDFLAGS  =  ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o -lpthread -lm -lbz2 -llzma -lz


all: libgab/libgab.a tabixpp/tabix.o samtools/bedidx.o bamtools/src/bamtools_fasta.o glactools 

%.o: %.cpp libgab/libgab.a tabixpp/tabix.o samtools/bedidx.o bamtools/src/bamtools_fasta.o 
	${CXX} ${CXXFLAGS} $< -o $@

libgab/gzstream/gzstream.o: libgab/libgab.a
	echo ""

libgab/libgab.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git

libgab/libgab.a: htslib/libhts.so libgab/libgab.h
	cd libgab/ && make libgab.a && make -C gzstream/  && cd ..

tabixpp/tabix.hpp:
	rm -rf tabixpp/
	git clone --recursive https://github.com/grenaud/tabixpp.git

tabixpp/tabix.o: tabixpp/tabix.hpp htslib/libhts.so
	make -C tabixpp

bamtools/src/bamtools_fasta.o:
	cd bamtools/src/ && g++ -lm -c -I.  -Ishared/ -Iutils/ utils/bamtools_fasta.cpp && cd ../..

samtools/bedidx.o: samtools/libbam.a
	echo ""

samtools/libst.a: samtools/libbam.a
	echo ""

samtools/libbam.a:  samtools/sam.h
	cd samtools/ && make && cd ..

samtools/sam.h: htslib/libhts.so
	rm -rf samtools/
	git clone --recursive https://github.com/samtools/samtools.git samtools/

htslib/libhts.a: htslib/libhts.so
	echo ""

htslib/libhts.so:  htslib/hts_internal.h
	cd htslib/ && git submodule update --init --recursive &&  make && cd ../

htslib/hts_internal.h:
	rm -rf htslib/
	git clone --recursive https://github.com/samtools/htslib.git

lib/libglactools.a: tabixpp/tabix.o libgab/libgab.a bamtools/src/bamtools_fasta.o  samtools/libbam.a samtools/libst.a GlacIndex.o GlacIDXSTATS.o  RandomGenomicCoord.o GenomicWindows.o GlacWindows.o GlacWriter.o MultiVCFreader.o GLF2ACF.o VcfMulti2ACF.o T3andme2ACF.o AXT2ACF.o Vcf2ACF.o Vcf2GLF.o  SingleAllele.o BAM2ACF.o EIGENSTRAT2ACF.o BEAGLE2GLF.o BPLINK2ACF.o AlleleRecords.o SingleGL.o GlacUndef.o  GlacMindef.o GlacSegsite.o GlacRename.o GlacSharing.o GlacNoSharing.o GlacNoStrictSharing.o GlacNosingle.o GlacBedfilter.o Glac2FREQSPEC.o Glac2BED.o Glac2VCF.o ACF2BPLINK.o ACF2FASTA.o ACF2GPHOCS.o ACF2NEXUS.o ACF2TREEMIX.o ACF2GROSS.o ACF2EIGENSTRAT.o ACF2BETASCAN.o GlacViewer.o GlacParser.o GlacStats.o GlacReplaceAncestor.o GlacUsePopAsRootAnc.o GlacCAT.o VCFreader.o SimpleVCF.o CoreVCF.o ReadTabix.o SetVCFFilters.o GlacMeld.o GlacDown.o GlacPopsub.o GlacCompute.o SumStatD.o DstatResult.o DstatCounter.o Dstat_core.o SumStatF3.o F3Result.o F3Counter.o F3_core.o SumStatF2.o F2Result.o F2Counter.o F2_core.o SumStatAvgCoa.o AvgCoaAlleleCounter.o AvgCoaResult.o ComputeAvgCoa_core.o ComputeFst_core.o DistAlleleCounter.o DistResult.o SumStatDist.o ComputeDist_core.o FstAlleleCounter.o FstResult.o SumStatFst.o GlacClosest.o GlacRemovepop.o GlacIntersect.o GlacUnion.o GlacReheader.o FilterVCF.o GlactoolsOperations.o GenomicRange.o  libgab/gzstream/gzstream.o 
	ar rs lib/libglactools.a $^

glactools:	glactools.o  lib/libglactools.a htslib/libhts.a samtools/bedidx.o 
	${CXX} -o $@ $^ $(LDLIBS)  $(LDFLAGS)

cleanall :
	make -C htslib/ clean
	make -C samtools/ clean
	make -C libgab/ clean
	make -C tabixpp/ clean
	make -C tabixpp/ clean
	rm bamtools/src/bamtools_fasta.o
	rm -f *.o glactools

clean :
	rm -f *.o glactools

.PHONY: test

test:	all
	cd test/ && bash testAll.sh && cd ..
