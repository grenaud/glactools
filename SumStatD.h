/*
 * SumStatD
 * Date: Oct-14-2016
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SumStatD_h
#define SumStatD_h


//#include "DstatCounter.h"
#include "DstatResult.h"
#include "Dstat_core.h"

/* #include "AvgCoaResult.h" */
/* #include "ComputeAvgCoa_core.h" */
#include "AlleleRecords.h"
#include "GlacParser.h"

//#include <signal.h>

using namespace std;

class SumStatD{
 private:
    DstatResult *** dstatResults;//[numberOfPopulations][numberOfPopulations][numberOfPopulations];
    unsigned int numberOfPopulations;
    bool onlyOnRef;
 public:
    vector<string> * populationNames ;
    DstatResult const * const * const * getDstatResult() const;

    SumStatD();
    SumStatD(const vector<string> * popNames,const bool onlyOnRef_);
    SumStatD(const SumStatD & other);
    ~SumStatD();
    SumStatD & operator=  (const SumStatD & other);


    SumStatD & operator+= (const SumStatD & other){
	if(other.numberOfPopulations != numberOfPopulations){
	    cerr<<"SumStatD.cpp  problem in =+ operator, adding different objects "<<other.numberOfPopulations<<" "<<numberOfPopulations<<endl;
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
		for(unsigned k=0;k<numberOfPopulations;k++){	       
		    dstatResults[i][j][k]+=other.dstatResults[i][j][k];
		}
	    }
	}

	return *this;
    }


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
		for(unsigned k=0;k<numberOfPopulations;k++){	       
		    dstatResults[i][j][k]-=other.dstatResults[i][j][k];
		}
	    }
	}
	
	return *this;
    }

    //void computeStat(    vector < AlleleRecords * >  * dataToUse,vector<string> & popNames,bool allowUndefined);
    void computeStat(       const   vector < AlleleRecords  >  * dataToUse,const vector<string> * popNames,bool allowUndefined);
    void computeStatSingle( const   AlleleRecords   * recordToUse                                         ,const bool allowUndefined);

    string print() const;
    //istream & read(istream & s);
    void read(const string & res);
    string printWithBootstraps(const   vector<SumStatD *> * bootstr,const string & dnaDistMode) const;

    friend ostream& operator<<(ostream& os, const SumStatD & ct){
    	os<<ct.print();
    	return os;
    }

    //TODO fix
    //istream& operator>>(istream& s,  SumStatDist& c){
    /* istream& operator>>(istream& s){ */
    /* 	cout<<"istream"<<endl; */
    /*    return read(s); */
    /* } */

};
#endif
