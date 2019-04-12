/*
 * SumStatPi
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatPi_h
#define SumStatPi_h
#include "ComputePi_core.h"
#include "AlleleRecords.h"
#include "GlacParser.h"
#include "PiResult.h"

using namespace std;

class SumStatPi{
 private:
    PiResult  ** pianceResults;//[numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;

 public:
    vector<string> * populationNames ;
    PiResult const * const * getPiResult() const;
    
    SumStatPi();
    SumStatPi(const vector<string> * popNames);
    SumStatPi(const SumStatPi & other);
    ~SumStatPi();
    SumStatPi & operator=  (const SumStatPi & other);

    SumStatPi & operator+= (const SumStatPi & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatPi.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatPi.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
	   
		piResults[i][j]+=other.piResults[i][j];

	    }
	}

	return *this;
    }

    /* SumStatPi operator-(const SumStatPi & s1){ */
    /* 	SumStatPi toReturn (s1); */
    /* 	//toReturn=-s2; */
    /* 	return toReturn; */
    /* } */

    SumStatPi & operator-= (const SumStatPi & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatPi.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatPi.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
		
		pianceResults[i][j]-=other.pianceResults[i][j];

	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat( const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined);

    string print() const;
    string printWithBootstraps(const   vector<SumStatPi *> * bootstr) const;

    friend ostream& operator<<(ostream& os, const SumStatPi & ct){
	os<<ct.print();
	return os;
    }

};
#endif
