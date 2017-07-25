/*
 * FilterVCF
 * Date: Aug-22-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FilterVCF.h"



static int rejectIndel=0;
static int rejectREFValidREF=0;
static int rejectREFValidALT=0;
static int rejectLOWCOV_REF=0;
static int rejectMap20=0;
static int rejectLOWMQ=0;
static int rejectLOWQUAL=0;
static int rejectCloseIndels=0;
static int rejectSysERR=0;
static int rejectRM=0;
static int rejectREF_unknownGeno=0;

//! To check whether a SimpleVCF passes the some filters
/*!
 *
 * This subroutine checks if a SimpleVCF record passes the following tests:
 * 
 *  0) The variation is not an indel
 *  1) Has A,C,G or T as reference allele
 *  2) Has A,C,G,T or . as alternative allele
 *  3) Has sequence coverage (depth) between two given values
 *  4) If the Map20 field is present, if checks if it is above a given cutoff
 *  5) Has the root mean square of the mapping quality score above a certain cutoff
 *  6) Has GQ field above a certain cutoff
 *  7) If not flagged as close to indel
 *  8) If not flagged as  syserr (systematic error)
 *  9) If not flagged as rm (repeat masked)\
 * 10) If the genotype is unknown (GT field = ./.)
 *
 \param smvcf               The SimpleVCF object to check, cannot be const because it can call the info field parser
 \param minCovcutoff        The minimum coverage cutoff
 \param maxCovcutoff        The maximum coverage cutoff
 \param minMapabilitycutoff    The mapability score cutoff
 \param minMQcutoff            The cutoff for the root mean square of the mapability score above a certain cutoff
 \param minGQcutoff            The cutoff for the GQ field (genotype quality)
 \return  : True if the SimpleVCF passes all the aforementioned checks
*/

bool passedFilters(SimpleVCF * smvcf,const SetVCFFilters * filtersToUse){
		   //int minCovcutoff,int maxCovcutoff,double minMapabilitycutoff,int minMQcutoff,int minGQcutoff){
    /////////////////////////////////////////////////////////////////
    //           FILTERING INDELS AND UNDEFINED REF ALLELE         //
    /////////////////////////////////////////////////////////////////
    
    if(!smvcf->isResolvedSingleBasePairREF()){ rejectREFValidREF++; 
	return false; }  //REF al  

    if(!smvcf->isResolvedSingleBasePairALT()){ 
	rejectREFValidALT++; 
	return false;  
    } //ALT al


    bool notFoundNonindel=smvcf->containsIndel() ;

    if(notFoundNonindel){ 
	rejectIndel++;  //we don't want indels		
	return false;
    }
    

    //rejecting unknown genotype (GT= ./.)
    if(smvcf->isUnresolvedGT() ){ 
	rejectREF_unknownGeno++; 
	return false; 
    }

    //if the SetVCFFilters says that we do not filter beyond that, return true
    if(filtersToUse->getDonotFilter()){
	return true;
    }

    if(filtersToUse->getDonotFilterButMQ()){//just need to check MQ field
	
	double minMQ=smvcf->getInfoField<double>("MQ");
		    
	if(minMQ < filtersToUse->getMinMQcutoff()){
	    rejectLOWMQ++;
	    return false;
	}else{
	    return true;
	}

    }



    /////////////////////////////////
    //FILTERING ON GENOTYPE QUALITY /
    /////////////////////////////////
		  


    // - are in the 2.5% tails of the coverage distribution
    int coverageREF=smvcf->getDepth();
    //cerr<<"coverageREF\t"<<coverageREF<<endl;
    if(coverageREF == -1){
	coverageREF=smvcf->getDepthInfo();
    }
    //cerr<<"coverageREF\t"<<coverageREF<<endl;
    if(coverageREF < filtersToUse->getMinCovcutoff() ||
       coverageREF > filtersToUse->getMaxCovcutoff() ){
	rejectLOWCOV_REF++;
	return false;
    }





    // - have GQ < X	
    //cerr<<"qual "<<smvcf->getGenotypeQual()<<endl;
    
    if( (smvcf->getGenotypeQual()!=-1)
	&&
	(smvcf->getGenotypeQual() < filtersToUse->getMinGQcutoff()) ){
	rejectLOWQUAL++;
	return false;
    }


		    

	    
    // - are ± Xbp of InDels
    if(filtersToUse->getFilterIndelProx()){
	bool isCloseToIndel=(smvcf->getCloseIndel() );	    
	if(isCloseToIndel){	
	    rejectCloseIndels++;
	    return false;
	}
    }

    // - are flagged as SysErr
    if(filtersToUse->getSystemError()){
	if(smvcf->hasInfoField("SysErr")){
	    rejectSysERR++;
	    return false;
	}
    }

    // - are flagged as RM (repeat masked)
    if(filtersToUse->getRepeatMasking()){
	if(smvcf->hasInfoField("RM")){
	    rejectRM++;
	    return false;
	}       
    }


		    
		    

    //Check the info fields at the end
    //Because the info fields are not parsed automatically
    

    // - have Map20 < 1
    if(smvcf->hasInfoField("Map20")){
	if( (smvcf->getInfoField<double>("Map20")) < filtersToUse->getMinMapabilitycutoff() ){
	    rejectMap20++;
	    return false;
	}	
    }





    // - have MQ< 30
    double minMQ=smvcf->getInfoField<double>("MQ");
    //cerr<<minMQ<<endl;
    if(minMQ < filtersToUse->getMinMQcutoff()){
	rejectLOWMQ++;
	return false;
    }

    return true;		   
}

string rejectFiltersTally(){

    return 
	"rejectIndel            "+stringify(rejectIndel)+"\n"+
	"rejectREFValidREF      "+stringify(rejectREFValidREF)+"\n"+
	"rejectREFValidALT      "+stringify(rejectREFValidALT)+"\n"+
	"rejectLOWCOV_REF       "+stringify(rejectLOWCOV_REF)+"\n"+
	"rejectMap20            "+stringify(rejectMap20)+"\n"+
	"rejectLOWMQ            "+stringify(rejectLOWMQ)+"\n"+
	"rejectLOWQUAL          "+stringify(rejectLOWQUAL)+"\n"+
	"rejectCloseIndels      "+stringify(rejectCloseIndels)+"\n"+
	"rejectSysERR           "+stringify(rejectSysERR)+"\n"+
	"rejectRM               "+stringify(rejectRM)+"\n"+
	"rejectREF_unknownGeno  "+stringify(rejectREF_unknownGeno)+"\n";


}
