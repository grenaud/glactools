/*
 * glactools
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "libgab.h"

#include "BAM2ACF.h"
#include "EIGENSTRAT2ACF.h"
#include "BEAGLE2GLF.h"
#include "BPLINK2ACF.h"

#include "Vcf2ACF.h"
#include "GLF2ACF.h"
#include "Vcf2GLF.h"
#include "VcfMulti2ACF.h"
#include "T3andme2ACF.h"
#include "AXT2ACF.h"

#include "GlacClosest.h"
#include "Glac2VCF.h"
#include "Glac2BED.h"
#include "GlacMeld.h"
#include "GlacDown.h"
#include "GlacPopsub.h"
#include "GlacRemovepop.h"
#include "GlacReheader.h"
#include "GlacCompute.h"
#include "GlacReplaceAncestor.h"
#include "GlacUsePopAsRootAnc.h"

#include "GlacIntersect.h"
#include "GlacUnion.h"
#include "Glac2FREQSPEC.h"


#include "GlacUndef.h"
#include "GlacMindef.h"

#include "GlacBedfilter.h"
#include "GlacSegsite.h"
#include "GlacRename.h"


#include "GlacSharing.h"
#include "GlacNoSharing.h"
#include "GlacNoStrictSharing.h"
#include "GlacNosingle.h"

#include "GlacCAT.h"
#include "ACF2BPLINK.h" 
#include "ACF2FASTA.h"
#include "ACF2GPHOCS.h"
#include "ACF2NEXUS.h"
#include "ACF2TREEMIX.h" 
#include "ACF2GROSS.h" 
#include "ACF2EIGENSTRAT.h" 
#include "ACF2BETASCAN.h" 

#include "GlacIndex.h"
#include "GlacIDXSTATS.h"
#include "GlacWindows.h"
#include "GlacStats.h"
#include "GlacViewer.h"

using namespace std;

int main (int argc, char *argv[]) {

    const string usage=string("\n")+
	"   This program contains various functions for managing:\n"+                       
	"   Genotype Likelihoods & Allele Counts\n"+                       
	"\n"+
	"   glf   genotype likelihood format\n"+
	"   acf   allele count format\n"+			      
	"\n"+                       
	"   "+string(argv[0]) +" <command> options\n"+                       
	"\n"+                       
	"   Commands:\n"+                       
	"   --Data import:\n"+                       
	"      vcf2acf         Convert single sample VCF to acf "+"\n"+
	"      vcf2glf         Convert single sample VCF to glf "+"\n"+
	"      vcfm2acf        Convert multi  sample VCF to acf "+"\n"+
	"      bam2acf         Convert single sample BAM to acf "+"\n"+
	"      axt2acf         Convert AXT alignment to acf "+"\n"+
	"      eigen2acf       Convert EIGENSTRAT data to acf "+"\n"+
	"      beagle2glf      Convert BEAGLE data to glf "+"\n"+
	"      bplink2acf      Convert binary PLINK to acf"+"\n"+
	"      23andme2acf     Convert 23andme data to acf "+"\n"+
	
	"\n"+
	"   --Filtering:\n"+                       
	"      noundef         No undefined sites for populations"+"\n"+
	"      mindef          Require all ind/pop to carry a min. # of defined alleles"+"\n"+
	"      bedfilter       Filter ACF/GLF file using sorted bedfile"+"\n"+
	"      segsite         Just retain segregating sites (or trans./transi)"+"\n"+
	"      sharing         Retain sites that share alleles between populations"+"\n"+
	"      nosharing       Retain sites that do NOT share alleles between populations"+"\n"+
	"      snosharing      Retain sites that STRICKLY do NOT share alleles between populations"+"\n"+
	"      nosingle        Retain variant sites where 2 individuals/pop or more have the variant"+"\n"+
	"\n"+
	"   --Computations:\n"+                       
	"      freqspec        Compute the frequency spectrum"+"\n"+
	"      closest         Return the distance between records"+"\n"+
	"      compute         Compute summary statistics"+"\n"+
	"      stats           Provide very basic stats"+"\n"+

	"\n"+
	"   --File transformations:\n"+                       
	"      cat             Concatenate (GL|AC)f files"+"\n"+
	"      intersect       Intersection of (GL|AC)f files"+"\n"+
	"      union           Union of (GL|AC)f files"+"\n"+
	"      reheader        Replace header of (GL|AC)f files"+"\n"+

	"\n"+
	"   --Population transformations:\n"+                       
	"      meld            Merge multiple populations as one for ACF files\n"+
	"      popsub          Keep a subset of the populations\n"+
	"      removepop       Remove a subset of the populations\n"+
	"      replaceanc      Use ancestral/root information from another file\n"+
	"      usepopsrootanc  Use 2 specified pops as ancestral/root information\n"+
	"      rename          Rename populations\n"+
	"      down            Downsample the allele count of all/some populations\n"+
	"\n"+
	"   --Data export:\n"+                       
	"      glac2bed        Convert a (GL|AC)f file to BED"+"\n"+
	"      glac2vcf        Convert a (GL|AC)f file to VCF"+"\n"+
	"      acf2bplink      Convert an ACF file to binary PLINK"+"\n"+
	"      acf2fasta       Convert an ACF file to fasta"+"\n"+ 
	"      acf2gphocs      Convert an ACF file to G-PhoCs"+"\n"+
	"      acf2nexus       Convert an ACF file to Nexus"+"\n"+
	"      acf2treemix     Convert an ACF file to treemix"+"\n"+
	"      acf2gross       Convert an ACF file to GRoSS"+"\n"+
	"      acf2eigenstrat  Convert an ACF file to EIGENSTRAT"+"\n"+
	"      acf2betascan    Convert an ACF file to betascan"+"\n"+
	//"      vcfm2glf      Convert multi  sample VCF to glf "+"\n"+
	"\n"+
	"   --GLF/ACF conversion:\n"+                       
	"      glf2acf         Convert glf to acf "+"\n"+
	"\n"+
	"   --Indexing:\n"+                       
	"      index           Index acf/glf file"+"\n"+
	"      idxstats        Basic statistics using the index of a acf/glf file"+"\n"+
	"\n"+
	"   --Viewing:\n"+                       
	"      view            View all or a region of a ACF/GLF file "+"\n"+
	"      windows         Print windows along the genome "+"\n"+

	"";

                              
    if( argc==1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "Usage "<<usage<<endl;
        return 1;       
    }

    if(string(argv[1]) == "vcf2acf"){
	Vcf2ACF  vcf2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<vcf2acf_.usage()<<endl;
	    return 1;       
	}
	argv++;
	argc--;
	return vcf2acf_.run(argc, argv);
	

    }else{      if(string(argv[1]) == "vcf2glf"){
	Vcf2GLF  vcf2glf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<vcf2glf_.usage()<<endl;
	    return 1;       
	}
	argv++;
	argc--;
	return vcf2glf_.run(argc, argv);

    }else{      if(string(argv[1]) == "vcfm2acf"){
	VcfMulti2ACF  vcfmulti2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<vcfmulti2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return vcfmulti2acf_.run(argc, argv);

    }else{      if(string(argv[1]) == "bam2acf"){
	BAM2ACF  bam2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<bam2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return bam2acf_.run(argc, argv);

    }else{      if(string(argv[1]) == "axt2acf"){
	AXT2ACF  axt2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<axt2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return axt2acf_.run(argc, argv);
		    
    }else{      if(string(argv[1]) == "eigen2acf"){
	EIGENSTRAT2ACF  eigen2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<eigen2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return eigen2acf_.run(argc, argv);
	
    }else{      if(string(argv[1]) == "beagle2glf"){
	BEAGLE2GLF  beagle2glf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<beagle2glf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return beagle2glf_.run(argc, argv);
	
    }else{      if(string(argv[1]) == "bplink2acf"){
	BPLINK2ACF  bplink2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<bplink2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return bplink2acf_.run(argc, argv);
		    
    }else{      if(string(argv[1]) == "glf2acf"){
	GLF2ACF  glf2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glf2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glf2acf_.run(argc, argv);
		    
	    	    
    }else{      if(string(argv[1]) == "23andme2acf"){
	T3andme2ACF  t3andme2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<t3andme2acf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return t3andme2acf_.run(argc, argv);

    }else{      if(string(argv[1]) == "glac2bed"){
	Glac2BED  glac2bed_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glac2bed_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glac2bed_.run(argc, argv);

    }else{      if(string(argv[1]) == "glac2vcf"){
	Glac2VCF  glac2vcf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glac2vcf_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glac2vcf_.run(argc, argv);

    }else{      if(string(argv[1]) == "closest"){
     	GlacClosest  glacclosest_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacclosest_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacclosest_.run(argc, argv);
		    	    	    
    }else{      if(string(argv[1]) == "acf2bplink"){
	ACF2BPLINK  acf2bplink_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2bplink_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2bplink_.run(argc, argv);

    }else{      if(string(argv[1]) == "acf2fasta"){
	ACF2FASTA  acf2fasta_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2fasta_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2fasta_.run(argc, argv);

    }else{      if(string(argv[1]) == "acf2gphocs"){
	ACF2GPHOCS  acf2gphocs_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2gphocs_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2gphocs_.run(argc, argv);


    }else{      if(string(argv[1]) == "acf2nexus"){
	ACF2NEXUS  acf2nexus_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2nexus_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2nexus_.run(argc, argv);

    }else{      if(string(argv[1]) == "acf2treemix"){
	ACF2TREEMIX  acf2treemix_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2treemix_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2treemix_.run(argc, argv);

    }else{      if(string(argv[1]) == "acf2betascan"){
	ACF2BETASCAN  acf2betascan_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2betascan_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2betascan_.run(argc, argv);

    }else{      if(string(argv[1]) == "acf2gross"){
	ACF2GROSS  acf2gross_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2gross_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2gross_.run(argc, argv);

    }else{      if(string(argv[1]) == "acf2eigenstrat"){
	ACF2EIGENSTRAT  acf2eigenstrat_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<acf2eigenstrat_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return acf2eigenstrat_.run(argc, argv);


    }else{      if(string(argv[1]) == "noundef"){
	GlacUndef  glacundef_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacundef_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacundef_.run(argc, argv);

    }else{      if(string(argv[1]) == "mindef"){
	GlacMindef  glacmindef_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacmindef_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacmindef_.run(argc, argv);

    }else{      if(string(argv[1]) == "bedfilter"){
	GlacBedfilter  glacbedfilter_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacbedfilter_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacbedfilter_.run(argc, argv);

    }else{      if(string(argv[1]) == "segsite"){
	GlacSegsite  glacsegsite_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacsegsite_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacsegsite_.run(argc, argv);

    }else{      if(string(argv[1]) == "rename"){
	GlacRename  glacrename_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacrename_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacrename_.run(argc, argv);

    }else{      if(string(argv[1]) == "sharing"){
	GlacSharing  glacsharing_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacsharing_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacsharing_.run(argc, argv);

    }else{      if(string(argv[1]) == "nosharing"){
	GlacNoSharing  glacnosharing_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacnosharing_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacnosharing_.run(argc, argv);
    }else{      if(string(argv[1]) == "snosharing"){
	GlacNoStrictSharing  glacnostrictsharing_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacnostrictsharing_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacnostrictsharing_.run(argc, argv);
	
    }else{      if(string(argv[1]) == "nosingle"){
	GlacNosingle  glacnosingle_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacnosingle_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacnosingle_.run(argc, argv);

    }else{      if(string(argv[1]) == "freqspec"){
	Glac2FREQSPEC  glac2freqspec_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glac2freqspec_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glac2freqspec_.run(argc, argv);

    }else{      if(string(argv[1]) == "cat"){
	GlacCAT  glaccat_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glaccat_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glaccat_.run(argc, argv);

    }else{      if(string(argv[1]) == "meld"){
	GlacMeld  glacmeld_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacmeld_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacmeld_.run(argc, argv);

    }else{      if(string(argv[1]) == "down"){
	GlacDown  glacdown_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacdown_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacdown_.run(argc, argv);

    }else{      if(string(argv[1]) == "popsub"){
	GlacPopsub  glacpopsub_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacpopsub_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacpopsub_.run(argc, argv);

    }else{      if(string(argv[1]) == "removepop"){
	GlacRemovepop  glacremovepop_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacremovepop_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacremovepop_.run(argc, argv);

    }else{      if(string(argv[1]) == "replaceanc"){
	GlacReplaceAncestor  glacreplaceanc_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacreplaceanc_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacreplaceanc_.run(argc, argv);

    }else{      if(string(argv[1]) == "usepopsrootanc"){
	GlacUsePopAsRootAnc  glacusepopsrootanc_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacusepopsrootanc_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacusepopsrootanc_.run(argc, argv);

    }else{      if(string(argv[1]) == "stats"){
	GlacStats  glacstats_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacstats_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacstats_.run(argc, argv);

    }else{      if(string(argv[1]) == "removepop"){
	GlacRemovepop  glacremovepop_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacremovepop_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacremovepop_.run(argc, argv);

    }else{      if(string(argv[1]) == "intersect"){
	GlacIntersect  glacintersect_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacintersect_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacintersect_.run(argc, argv);


    }else{      if(string(argv[1]) == "union"){
	GlacUnion  glacunion_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacunion_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacunion_.run(argc, argv);

    }else{      if(string(argv[1]) == "reheader"){
	GlacReheader  glacreheader_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacreheader_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacreheader_.run(argc, argv);

    }else{      if(string(argv[1]) == "compute"){
	GlacCompute  glaccompute_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glaccompute_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glaccompute_.run(argc, argv);

		    	    	    
    }else{      if(string(argv[1]) == "index"){
	GlacIndex  glindex_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glindex_.usage()<<endl;
	    return 1;       
	}


	argv++;
	argc--;
	return glindex_.run(argc, argv);
		    
	    

    }else{      if(string(argv[1]) == "idxstats"){
	GlacIDXSTATS  glacIdxstats_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacIdxstats_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacIdxstats_.run(argc, argv);		   	    

    }else{      if(string(argv[1]) == "windows"){
	GlacWindows  glacWindows_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glacWindows_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glacWindows_.run(argc, argv);		   	    

    }else{      if(string(argv[1]) == "view"){
	GlacViewer  glviewer_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<glviewer_.usage()<<endl;
	    return 1;       
	}

	argv++;
	argc--;
	return glviewer_.run(argc, argv);
		    
	    
    }else{      	    
	cerr<<"invalid command "<<string(argv[1])<<endl;
	return 1;
																								}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
    
    return 0;
}

