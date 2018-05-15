/*
 * DistAlleleCounter
 * Date: Aug-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef DistAlleleCounter_h
#define DistAlleleCounter_h

#include <iostream>
#include <ostream> 
#include <math.h>
#include <limits>
//#include <iostream>
//#include <fstream>
using namespace std;

class DistAlleleCounter{
 private:

 public:
    void reinitializedCounters();

    //counter of mutations
    unsigned int counterSame;
    unsigned int counterCommon;
    unsigned int counterReference;
    unsigned int counterSample;

    DistAlleleCounter(); //constructor

    double  lowerConf  (const unsigned int shortBranch,const unsigned int commonBranch) const;
    double  highConf   (const unsigned int shortBranch,const unsigned int commonBranch) const;
    /* friend ostream & operator<<(ostream & os, const DistAlleleCounter & ct); */
    string getHeader(string prefixToAdd="");
    double distRefSam () const ;
    double distSamRef () const ;

    DistAlleleCounter & operator+= (const DistAlleleCounter & other){
	counterSame      += other.counterSame;
	counterCommon    += other.counterCommon;
	counterReference += other.counterReference;
	counterSample    += other.counterSample;
	return *this;
    }

    DistAlleleCounter & operator-= (const DistAlleleCounter & other){
	counterSame      -= other.counterSame;
	counterCommon    -= other.counterCommon;
	counterReference -= other.counterReference;
	counterSample    -= other.counterSample;
	return *this;
    }


    friend ostream& operator<<(ostream& os, const DistAlleleCounter & ct){
	os<<ct.counterSame<<"\t"
	  <<ct.counterCommon<<"\t"
	  <<ct.counterReference<<"\t"
	  <<ct.counterSample<<"\t";
	if( (ct.counterReference+ct.counterCommon) != 0)	
	    os<<float(ct.counterReference)/float(ct.counterReference+ct.counterCommon)<<"\t"<<ct.lowerConf(ct.counterReference,ct.counterCommon)<<"\t"<<ct.highConf(ct.counterReference,ct.counterCommon);       
	else
	    os<<"nan"<<"\t"<<"nan"<<"\t"<<"nan";
	os<<"\t";
	if( (ct.counterSample+ct.counterCommon) != 0)
	    os<<float(ct.counterSample)/float(ct.counterSample+ct.counterCommon)<<"\t"<<ct.lowerConf(ct.counterSample,ct.counterCommon)<<"\t"<<ct.highConf(ct.counterSample,ct.counterCommon);       
	else
	    os<<"nan"<<"\t"<<"nan"<<"\t"<<"nan";
	
	return os;
    }

};



#endif
