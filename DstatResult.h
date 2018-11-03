/*
 * DstatResult
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef DstatResult_h
#define DstatResult_h

using namespace std;

#include <list>
#include <vector>

#include "DstatCounter.h"

class DstatResult{
 private:
 public:

    //This is everything
    DstatCounter all;
    //This is without the ones marked as CpG
    DstatCounter noCpg;
    //This is only with the ones marked as CpG
    DstatCounter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    DstatCounter transversions;
    DstatCounter transitions;
    //This excludes the following cases:
    // S = T, A = C
    // S = A, A = G
    DstatCounter noDamage;

    DstatResult();
    DstatResult(const DstatResult & other);
    ~DstatResult();
    /* DstatResult & operator= (const DstatResult & other); */

    DstatResult & operator= (const DstatResult & other){
	this->all             = other.all;
	this->noCpg           = other.noCpg;
	this->onlyCpg         = other.onlyCpg;
	this->transversions   = other.transversions;
	this->transitions     = other.transitions;
	this->noDamage        = other.noDamage;
	return *this;
    }
    void addIfNotInfinity( vector<double> * target , double val );
    string printWithBootstrap(list<vector< vector< vector<DstatResult> >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int k,unsigned int numberOfBootstraps);
    string printWithJacknife(const vector< const DstatResult * >  * jacknife);


    DstatResult &  operator+=(const DstatResult & other);
    DstatResult &  operator-=(const DstatResult & other);


    friend ostream& operator<<(ostream& os, const DstatResult & dr){
	/* cout<<"DstatResult print() begin"<<endl; */
	/* exit(1); */
	os<<dr.all          <<"\t"
	  <<dr.onlyCpg      <<"\t"  
	  <<dr.noCpg        <<"\t"  
	  <<dr.transitions  <<"\t"  
	  <<dr.transversions<<"\t"  
	  <<dr.noDamage;  

	/* os<<":\t"               <<dr.all.headerForCount()  <<"\n"; */
	/* os<<"all_sites:\t"      <<dr.all                   <<"\n"; */
	/* os<<"nocpg_sites:\t"    <<dr.noCpg                 <<"\n"; */
	/* os<<"onlycpg_sites:\t"  <<dr.onlyCpg               <<"\n"; */
	/* os<<"transitions:\t"    <<dr.transitions           <<"\n"; */
	/* os<<"transversions:\t"  <<dr.transversions         <<"\n"; */
	/* os<<"nodamage:\t"       <<dr.noDamage              <<"\n"; */

	
	/* os<<"all\t"<<dr.all<<"\n" */
	/*   <<"all\t"<<dr.onlyCpg<<"\n" */
	/*   <<dr.noCpg<<"\n" */
	/*   <<dr.transitions<<"\n" */
	/*   <<dr.transversions<<"\n" */
	/*   <<dr.noDamage; */
	return os;
    }

    friend istream &operator>>(istream &is , DstatResult &dr){
	//cerr<<" op TEST"<<endl;
	is
	    >>dr.all
	    >>dr.onlyCpg
	    >>dr.noCpg
	    >>dr.transitions
	    >>dr.transversions
	    >>dr.noDamage;
	return is;
	//return dsr.read(s);
    }

    /* istream read(istream &s){ */
    /* 	cout<< */
};//class DstatResult;


///}//class DstatResult;


#endif
