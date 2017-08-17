/*
 * SetVCFFilters
 * Date: Jan-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 *
 *  This class describes a set of VCF filters
 * 
 */

#ifndef SetVCFFilters_h
#define SetVCFFilters_h

#include "utils.h"

using namespace std;

class SetVCFFilters{
private:
    int minGQcutoff;
    int minMQcutoff;
    double minMapabilitycutoff;
    bool filterIndelProx;
    bool repeatMasking  ;
    bool systemError  ;
    int  minCovcutoff  ;
    int  maxCovcutoff  ;
    string name;
    bool donotFilter  ;
    bool donotFilterButMQ  ;

public:
   
    SetVCFFilters(int minGQcutoff=40,
		  int minMQcutoff=30,
		  double minMapabilitycutoff=1.0,
		  bool filterIndelProx = true,
		  bool repeatMasking   = true,   
		  bool systemError     = true,
		  int  minCovcutoff  =0,
		  int  maxCovcutoff  =1000,
		  bool   donotFilter=false,
		  bool   donotFilterButMQ=false);
    ~SetVCFFilters();

    void setName(string name) ;
    string getName() const;

    int getMinGQcutoff() const ;
    int getMinMQcutoff() const ;
    double getMinMapabilitycutoff() const ;
    bool getFilterIndelProx() const ;
    bool getRepeatMasking() const ;
    bool getSystemError() const ;
    int  getMinCovcutoff() const ;
    int  getMaxCovcutoff() const ;
    bool getDonotFilter() const ;
    bool getDonotFilterButMQ() const ;

    friend ostream& operator<<(ostream& os, const SetVCFFilters & ct){
	if(ct.name != "")
	    os<<"Cutoffs for "<<ct.name<< " : "<<endl;
	else
	    os<<"Cutoffs  : "<<endl;

	if(ct.donotFilterButMQ){
	    os<<"None, but MQ"<<endl;
	    os<<"Minimum RMS of the mapping qualities (MQ)               = "<<ct.minMQcutoff<<endl;
	    os<<"--------------"<<endl;
	}else{ 	
	    if(ct.donotFilter){
		os<<"None, all sites are allowed to proceed"<<endl;
		os<<"--------------"<<endl;
	    }else{
		os<<"Minimum genotype quality (GQ)                           = "<<ct.minGQcutoff<<endl;
		os<<"Minimum RMS of the mapping qualities (MQ)               = "<<ct.minMQcutoff<<endl;
		os<<"Minimum mapability                                      = "<<ct.minMapabilitycutoff<<endl;
		os<<"Filter sites close to an indel                          = "<<booleanAsString(ct.filterIndelProx)<<endl;
		/* if(ct.filterIndelProx){ */
		/*     os<<"Proximity (in bp) for an indel "<<ct.bpForIndels<<endl; */
		/* } */
		os<<"Filter sites close labeled as repeat masked (RM)        = "<<booleanAsString(ct.repeatMasking)<<endl;
		os<<"Filter sites close labeled as systematic error (SysErr) = "<<booleanAsString(ct.systemError)<<endl;
		
		os<<"Minimum coverage                                        = "<<ct.minCovcutoff<<endl;
		os<<"Minimum coverage                                        = "<<ct.maxCovcutoff<<endl;
		os<<"--------------"<<endl;
	    }
	}
	return os;
    }
	    

};
#endif
