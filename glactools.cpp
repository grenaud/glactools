/*
 * glactools
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"

#include "BAM2ACF.h"
#include "Vcf2ACF.h"
#include "GLF2ACF.h"
#include "Vcf2GLF.h"
#include "VcfMulti2ACF.h"
#include "T3andme2ACF.h"

#include "Glac2BED.h"
#include "GlacCAT.h"
#include "ACF2BPLINK.h" 
#include "ACF2FASTA.h"
#include "ACF2GPHOCS.h"
#include "ACF2NEXUS.h"
#include "ACF2TREEMIX.h" 

#include "GlacIndex.h"
#include "GlacViewer.h"

using namespace std;

int main (int argc, char *argv[]) {

    const string usage=string("\n")+
	"\tThis program contains various functions for managing:\n"+                       
	"\tGenotype Likelihoods & Allele Counts\n"+                       
	"\n"+
	"\tglf\tgenotype likelihood format\n"+
	"\tacf\tallele count format\n"+			      
	"\n"+                       
	"\t"+string(argv[0]) +" <command> options\n"+                       
	"\n"+                       
	"\tCommands:\n"+                       
	"\t--Data import:\n"+                       
	"\t\tvcf2acf\t\t\tConvert single sample VCF to acf "+"\n"+
	"\t\tvcf2glf\t\t\tConvert single sample VCF to glf "+"\n"+
	"\t\tvcfm2acf\t\tConvert multi  sample VCF to acf "+"\n"+
	"\t\tbam2acf\t\t\tConvert single sample BAM to acf "+"\n"+
	"\t\t23andme2acf\t\tConvert 23andme data to acf "+"\n"+
	"\n"+
	"\t--Data export:\n"+                       
	"\t\tglac2bed\t\tConvert a (GL|AC)f file to BED "+"\n"+
	"\t\tacf2binaryplink\t\tConvert an ACF file to binary PLINK "+"\n"+
	"\t\tacf2fasta\t\tConvert an ACF file to fasta "+"\n"+ 
	"\t\tacf2gphocs\t\tConvert an ACF file to G-PhoCs "+"\n"+
	"\t\tacf2nexus\t\tConvert an ACF file to Nexus "+"\n"+
	"\t\tacf2treemix\t\tConvert an ACF file to treemix "+"\n"+
	"\n"+
	"\t--File transformations:\n"+                       
	"\t\tglaccat\t\t\tConcatenate (GL|AC)f files"+"\n"+
	"\t\tglacintersect\t\tIntersection of (GL|AC)f files"+"\n"+
	"\t\tglacunion\t\tUnion of (GL|AC)f files"+"\n"+
	"\n"+
	"\t--Population transformations:\n"+                       
	"\t\tglacmeld\t\t"+"\n"+
	//"\t\tvcfm2glf\t\tConvert multi  sample VCF to glf "+"\n"+
	"\n"+
	"\t--GLF/ACF conversion:\n"+                       
	"\t\tglf2acf\t\t\tConvert glf to acf "+"\n"+
	"\n"+
	"\t--Indexing:\n"+                       
	"\t\tindex\t\t\tindex acf/glf file"+"\n"+
	"\n"+
	"\t--Viewing:\n"+                       
	"\t\tview\t\t\tview all or a region of a acf/glf file "+"\n"+

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

    }else{      if(string(argv[1]) == "glaccat"){
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
							    }}}}}}}}}}}}}}}
    
    return 0;
}

