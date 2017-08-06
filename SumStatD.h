/*
 * SumStatD
 * Date: Oct-14-2016
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatD_h
#define SumStatD_h

#include "ComputeAvgCoa_core.h"
#include "AlleleRecords.h"
#include "GlacParser.h"
#include "AvgCoaResult.h"

using namespace std;

class SumStatD{
 private:
    AvgCoaResult  ** divergenceResults;//[numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;

 public:
    vector<string> * populationNames ;
    AvgCoaResult const * const * getAvgCoaResult() const;
    
    SumStatD();
    SumStatD(const vector<string> * popNames);
    SumStatD(const SumStatD & other);
    ~SumStatD();
    SumStatD & operator=  (const SumStatD & other);

    SumStatD & operator+= (const SumStatD & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatD.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatD.cpp  problem in =+ operator, adding different objects"<<endl;
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

    /* SumStatD operator-(const SumStatD & s1){ */
    /* 	SumStatD toReturn (s1); */
    /* 	//toReturn=-s2; */
    /* 	return toReturn; */
    /* } */

    SumStatD & operator-= (const SumStatD & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatD.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatD.cpp  problem in =+ operator, adding different objects"<<endl;
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
    string printWithBootstraps(const   vector<SumStatD *> * bootstr) const;

    friend ostream& operator<<(ostream& os, const SumStatD & ct){
	os<<ct.print();
	return os;
    }

};
#endif
