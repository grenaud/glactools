/*
 * SumStatAvgCoa
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatAvgCoa_h
#define SumStatAvgCoa_h
#include "ComputeAvgCoa_core.h"
#include "AlleleRecords.h"
#include "GlacParser.h"
#include "AvgCoaResult.h"

using namespace std;

class SumStatAvgCoa{
 private:
    AvgCoaResult  ** divergenceResults;//[numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;

 public:
    vector<string> * populationNames ;
    AvgCoaResult const * const * getAvgCoaResult() const;
    
    SumStatAvgCoa();
    SumStatAvgCoa(const vector<string> * popNames);
    SumStatAvgCoa(const SumStatAvgCoa & other);
    ~SumStatAvgCoa();
    SumStatAvgCoa & operator=  (const SumStatAvgCoa & other);

    SumStatAvgCoa & operator+= (const SumStatAvgCoa & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatAvgCoa.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatAvgCoa.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
	   
		divergenceResults[i][j]+=other.divergenceResults[i][j];

	    }
	}

	return *this;
    }

    /* SumStatAvgCoa operator-(const SumStatAvgCoa & s1){ */
    /* 	SumStatAvgCoa toReturn (s1); */
    /* 	//toReturn=-s2; */
    /* 	return toReturn; */
    /* } */

    SumStatAvgCoa & operator-= (const SumStatAvgCoa & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatAvgCoa.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatAvgCoa.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
		
		divergenceResults[i][j]-=other.divergenceResults[i][j];

	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat( const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined);

    string print() const;
    string printWithBootstraps(const   vector<SumStatAvgCoa *> * bootstr) const;

    friend ostream& operator<<(ostream& os, const SumStatAvgCoa & ct){
	os<<ct.print();
	return os;
    }

};
#endif
