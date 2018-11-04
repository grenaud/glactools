/*
 * SumStatF2
 * Date: Oct-14-2016
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatF2_h
#define SumStatF2_h


//#include "DstatCounter.h"
#include "F2Result.h"
#include "F2_core.h"

/* #include "AvgCoaResult.h" */
/* #include "ComputeAvgCoa_core.h" */
#include "AlleleRecords.h"
#include "GlacParser.h"

//#include <signal.h>

using namespace std;

class SumStatF2{
 private:
    F2Result ** f2Results;//[numberOfPopulations][numberOfPopulations]
    unsigned int numberOfPopulations;
    //TODO does not need to be initialized at each site
    /* double average_nInds[ numberOfPopulations]; */
    /* double average_hzy[   numberOfPopulations]; */
    /* double id2nsnp[       numberOfPopulations]; */

    /* vector<double> average_nInds; */
    /* vector<double> average_hzy; */
    /* vector<double> id2nsnp; */

 public:
    vector<string> * populationNames ;
    F2Result const * const * getF2Result() const;

    SumStatF2();
    SumStatF2(const vector<string> * popNames);
    SumStatF2(const SumStatF2 & other);
    ~SumStatF2();
    SumStatF2 & operator=  (const SumStatF2 & other);


    SumStatF2 & operator+= (const SumStatF2 & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatF2.cpp  problem in =+ operator, adding different objects "<<other.numberOfPopulations<<" "<<numberOfPopulations<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatF2.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       	       
		f2Results[i][j] +=other.f2Results[i][j];
	    }
	}
	

	return *this;
    }


    SumStatF2 & operator-= (const SumStatF2 & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatF2.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatF2.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       		
		f2Results[i][j] -= other.f2Results[i][j];
	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat(       const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,      bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse                                         ,const bool allowUndefined);
    void read(const string & res);
    string print() const;
    string printWithBootstraps(const   vector<SumStatF2 *> * bootstr,const string & dnaDistMode) const;

    friend ostream& operator<<(ostream& os, const SumStatF2 & ct){
	os<<ct.print();
	return os;
    }

};
#endif
