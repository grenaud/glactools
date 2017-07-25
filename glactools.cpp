/*
 * glactools
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "vcf2acf.h"
#include "vcf2glf.h"

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
	"\t\tvcf2acf\t\tConvert VCF to acf "+"\n"+
	"\t\tvcf2glf\t\tConvert VCF to glf "+"\n"+
			      "";

                              
    if( argc==1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "Usage "<<usage<<endl;
        return 1;       
    }

    if(string(argv[1]) == "vcf2acf"){
	vcf2acf  vcf2acf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<vcf2acf_.usage()<<endl;
	    return 1;       
	}
	// char ** tosend =argv;

	// cout<<tosend[0]<<endl;
	// tosend++;
	// cout<<tosend[0]<<endl;
	argv++;
	argc--;
	return vcf2acf_.run(argc, argv);
	

    }else{      if(string(argv[1]) == "vcf2glf"){
	vcf2glf  vcf2glf_;

	if( argc==2 ||
	    (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
	    ){	    
	    cerr<<vcf2glf_.usage()<<endl;
	    return 1;       
	}
	// char ** tosend =argv;

	// cout<<tosend[0]<<endl;
	// tosend++;
	// cout<<tosend[0]<<endl;
	argv++;
	argc--;
	return vcf2glf_.run(argc, argv);
		    
	    
    }else{      
	    
	    cerr<<"invalid command "<<string(argv[1])<<endl;
	    return 1;
	}}
    
    return 0;
}

