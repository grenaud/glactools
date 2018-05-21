/*
 * SumStatDist
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatDist_h
#define SumStatDist_h
#include "ComputeDist_core.h"
#include "AlleleRecords.h"
#include "GlacParser.h"
#include "DistResult.h"

using namespace std;

class SumStatDist{
 private:
    DistResult  ** distanceResults;//[numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;

 public:
    vector<string> * populationNames ;
    DistResult const * const * getDistResult() const;
    
    SumStatDist();
    SumStatDist(const vector<string> * popNames);
    SumStatDist(const SumStatDist & other);
    ~SumStatDist();
    SumStatDist & operator=  (const SumStatDist & other);

    SumStatDist & operator+= (const SumStatDist & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatDist.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatDist.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
	   
		distanceResults[i][j]+=other.distanceResults[i][j];

	    }
	}

	return *this;
    }

    /* SumStatDist operator-(const SumStatDist & s1){ */
    /* 	SumStatDist toReturn (s1); */
    /* 	//toReturn=-s2; */
    /* 	return toReturn; */
    /* } */

    SumStatDist & operator-= (const SumStatDist & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatDist.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatDist.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
		
		distanceResults[i][j]-=other.distanceResults[i][j];

	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat( const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined);

    string print() const;
    string printWithBootstraps(const   vector<SumStatDist *> * bootstr,const string & dnaDistMode) const;

    friend ostream& operator<<(ostream& os, const SumStatDist & ct){
	os<<ct.print();
	return os;
    }

};
#endif
