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
			      "\t\tvcf2acf\t\tConvert VCF to either glf or acf "+
			      "";
                              
    if( argc==1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "Usage "<<usage<<endl;
        return 1;       
    }

    if(string(argv[1]) == "vcf2acf"){
	vcf2acf  vcf2acf_;
	cout<<vcf2acf_.usage()<<endl;
	

    }else{ 
	cerr<<"invalid command "<<string(argv[1])<<endl;
	return 1;
    }
    return 0;
}

