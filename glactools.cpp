/*
 * glactools
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "Vcf2ACF.h"
#include "GLF2ACF.h"
#include "Vcf2GLF.h"
#include "VcfMulti2ACF.h"

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
	"\t\tvcf2acf\t\t Convert single sample VCF to acf "+"\n"+
	"\t\tvcf2glf\t\t Convert single sample VCF to glf "+"\n"+
	"\t\tvcfm2acf\t\tConvert multi  sample VCF to acf "+"\n"+
	//"\t\tvcfm2glf\t\tConvert multi  sample VCF to glf "+"\n"+
	"\n"+
	"\t--GLF/ACF conversion:\n"+                       
	"\t\tglf2acf\t\tConvert glf to acf "+"\n"+
	"\n"+
	"\t--Indexing:\n"+                       
	"\t\tindex\t\tindex acf/glf file"+"\n"+
	"\n"+
	"\t--Viewing:\n"+                       
	"\t\tview\t\tview all or a region of a acf/glf file "+"\n"+

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
			}}}}}}
    
    return 0;
}

