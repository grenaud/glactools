/*
 * F3Result
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef F3Result_h
#define F3Result_h

using namespace std;

#include <list>
#include <vector>

#include "F3Counter.h"

class F3Result{
 private:
 public:

    //This is everything
    F3Counter all;
    //This is without the ones marked as CpG
    F3Counter noCpg;
    //This is only with the ones marked as CpG
    F3Counter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    F3Counter transversions;
    F3Counter transitions;
    //This excludes the following cases:
    // S = T, A = C
    // S = A, A = G
    F3Counter noDamage;

    F3Result();
    F3Result(const F3Result & other);
    ~F3Result();
    /* F3Result & operator= (const F3Result & other); */

    F3Result & operator= (const F3Result & other){
	this->all             = other.all;
	this->noCpg           = other.noCpg;
	this->onlyCpg         = other.onlyCpg;
	this->transversions   = other.transversions;
	this->transitions     = other.transitions;
	this->noDamage        = other.noDamage;
	return *this;
    }
    void addIfNotInfinity( vector<double> * target , double val );
    string printWithBootstrap(list<vector< vector< vector<F3Result> >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int k,unsigned int numberOfBootstraps);
    string printWithJacknife(const vector< const F3Result * >  * jacknife);


    F3Result &  operator+=(const F3Result & other);
    F3Result &  operator-=(const F3Result & other);

    friend ostream& operator<<(ostream& os, const F3Result & dr){
	/* cout<<"F3Result print() begin"<<endl; */
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

};
#endif
