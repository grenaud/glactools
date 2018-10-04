/*
 * SumStatF3
 * Date: Oct-14-2016
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatF3_h
#define SumStatF3_h


//#include "DstatCounter.h"
#include "F3Result.h"
#include "F3_core.h"

/* #include "AvgCoaResult.h" */
/* #include "ComputeAvgCoa_core.h" */
#include "AlleleRecords.h"
#include "GlacParser.h"

//#include <signal.h>

using namespace std;

class SumStatF3{
 private:
    F3Result *** f3Results;//[numberOfPopulations][numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;

 public:
    vector<string> * populationNames ;
    F3Result const * const * const * getF3Result() const;

    SumStatF3();
    SumStatF3(const vector<string> * popNames);
    SumStatF3(const SumStatF3 & other);
    ~SumStatF3();
    SumStatF3 & operator=  (const SumStatF3 & other);


    SumStatF3 & operator+= (const SumStatF3 & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatF3.cpp  problem in =+ operator, adding different objects "<<other.numberOfPopulations<<" "<<numberOfPopulations<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatF3.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
		for(unsigned k=0;k<numberOfPopulations;k++){	       
		    f3Results[i][j][k]+=other.f3Results[i][j][k];
		}
	    }
	}

	return *this;
    }


    SumStatF3 & operator-= (const SumStatF3 & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatF3.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatF3.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
		for(unsigned k=0;k<numberOfPopulations;k++){	       
		    f3Results[i][j][k]-=other.f3Results[i][j][k];
		}
	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat(       const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse                                         ,const bool allowUndefined);

    string print() const;
    string printWithBootstraps(const   vector<SumStatF3 *> * bootstr,const string & dnaDistMode) const;

    friend ostream& operator<<(ostream& os, const SumStatF3 & ct){
	os<<ct.print();
	return os;
    }

};
#endif
