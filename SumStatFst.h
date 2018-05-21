/*
 * SumStatFst
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatFst_h
#define SumStatFst_h
#include "ComputeFst_core.h"
#include "AlleleRecords.h"
#include "GlacParser.h"
#include "FstResult.h"

using namespace std;

class SumStatFst{
 private:
    FstResult  ** fstResults;//[numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;

 public:
    vector<string> * populationNames ;
    FstResult const * const * getFstResult() const;
    
    SumStatFst();
    SumStatFst(const vector<string> * popNames);
    SumStatFst(const SumStatFst & other);
    ~SumStatFst();
    SumStatFst & operator=  (const SumStatFst & other);

    SumStatFst & operator+= (const SumStatFst & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatFst.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatFst.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
	   
		fstResults[i][j]+=other.fstResults[i][j];

	    }
	}

	return *this;
    }

    /* SumStatFst operator-(const SumStatFst & s1){ */
    /* 	SumStatFst toReturn (s1); */
    /* 	//toReturn=-s2; */
    /* 	return toReturn; */
    /* } */

    SumStatFst & operator-= (const SumStatFst & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatFst.cpp  problem in =+ operator, adding different objects"<<endl;
	    exit(1);
	}

	for(unsigned int i=0;i<populationNames->size();i++){ 

	    if(other.populationNames->at(i) != populationNames->at(i)){
		cerr<<"SumStatFst.cpp  problem in =+ operator, adding different objects"<<endl;
		exit(1);		
	    }

	}


	for(unsigned i=0;i<numberOfPopulations;i++){
	    for(unsigned j=0;j<numberOfPopulations;j++){	       
		
		fstResults[i][j]-=other.fstResults[i][j];

	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat( const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined);

    string print() const;
    string printWithBootstraps(const   vector<SumStatFst *> * bootstr, const string & dnaDistMode) const;

    friend ostream& operator<<(ostream& os, const SumStatFst & ct){
	os<<ct.print();
	return os;
    }

};
#endif
