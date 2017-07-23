/*
 * vcf2acf
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "vcf2acf.h"



vcf2acf::vcf2acf(){

}

vcf2acf::~vcf2acf(){

}

string vcf2acf::usage() const{

    return string(string("") +"vcf2acf <options> [vcf file] [name sample] "+"\n"+
		  "\nThis program convert VCF files into acf (prints to the stdout)\n"+                       
		  "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
		  "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
		  "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
		  "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
		  "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
		  
		  "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
		  // "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 
		  
		  "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
		  "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
		  "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
		  "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
		  "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n");
    
    //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+


}


