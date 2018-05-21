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
#include <vector>
#include <sstream>

//#include <iostream>
//#include <fstream>
using namespace std;

class DistAlleleCounter{
 private:

 public:
    unsigned int count[16];
    //vector<string> dimers ;
    static const string dimers[16];

    void reinitializedCounters();

    void addAllelePair(int indexDimer);
    string headerForCount() const;
    unsigned int getIdent() const;
    unsigned int getMutations() const;

    DistAlleleCounter(); //constructor

    /* double  lowerConf  (const unsigned int shortBranch,const unsigned int commonBranch) const; */
    /* double  highConf   (const unsigned int shortBranch,const unsigned int commonBranch) const; */
    /* friend ostream & operator<<(ostream & os, const DistAlleleCounter & ct); */
    string getHeader(string prefixToAdd="");
    /* double distRefSam () const ; */
    /* double distSamRef () const ; */
    double dist(const string & dnaDistMode) const;

    DistAlleleCounter & operator+= (const DistAlleleCounter & other){
	for(int i=0;i<16;i++)
	    count[i] += other.count[i] ;
	/* counterSame      += other.counterSame; */
	/* counterCommon    += other.counterCommon; */
	/* counterReference += other.counterReference; */
	/* counterSample    += other.counterSample; */
	return *this;
    }

    DistAlleleCounter & operator-= (const DistAlleleCounter & other){
	for(int i=0;i<16;i++)
	    count[i] -= other.count[i] ;
	/* counterSame      -= other.counterSame; */
	/* counterCommon    -= other.counterCommon; */
	/* counterReference -= other.counterReference; */
	/* counterSample    -= other.counterSample; */
	return *this;
    }

    string printWithLabel() const;

    friend ostream& operator<<(ostream& os, const DistAlleleCounter & ct){
	/* for(int i=0;i<15;i++) */
	/*     os<<ct.counter[i]<<"\t"; */
	/* os<<ct.counter[15]<<"\t"; */


	unsigned int mutations=0;
	unsigned int ident    =0;

	for(int i=0;i<16;i++){
	    if( i == 0 || i == 5 || i == 10 || i == 15 ){
		ident+=     ct.count[i] ;
	    }else{
		mutations+= ct.count[i];	
	    }
	}

	for(int i=0;i<15;i++){
	    os<<ct.count[i]<<"\t";
	}
	os<<ct.count[15]<<"\t"<<ident<<"\t"<<mutations<<"\t"<<double(mutations) /double(mutations+ident);

	/* os<<ct.counterSame<<"\t" */
	/*   <<ct.counterCommon<<"\t" */
	/*   <<ct.counterReference<<"\t" */
	/*   <<ct.counterSample<<"\t"; */
	/* if( (ct.counterReference+ct.counterCommon) != 0)	 */
	/*     os<<float(ct.counterReference)/float(ct.counterReference+ct.counterCommon)<<"\t"<<ct.lowerConf(ct.counterReference,ct.counterCommon)<<"\t"<<ct.highConf(ct.counterReference,ct.counterCommon);        */
	/* else */
	/*     os<<"nan"<<"\t"<<"nan"<<"\t"<<"nan"; */
	/* os<<"\t"; */
	/* if( (ct.counterSample+ct.counterCommon) != 0) */
	/*     os<<float(ct.counterSample)/float(ct.counterSample+ct.counterCommon)<<"\t"<<ct.lowerConf(ct.counterSample,ct.counterCommon)<<"\t"<<ct.highConf(ct.counterSample,ct.counterCommon);        */
	/* else */
	/*     os<<"nan"<<"\t"<<"nan"<<"\t"<<"nan"; */
	
	return os;
    }

};


#endif
