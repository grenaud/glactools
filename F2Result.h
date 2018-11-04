/*
 * F2Result
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef F2Result_h
#define F2Result_h

using namespace std;

#include <list>
#include <vector>

#include "F2Counter.h"

class F2Result{
 private:
 public:

    //This is everything
    F2Counter all;
    //This is without the ones marked as CpG
    F2Counter noCpg;
    //This is only with the ones marked as CpG
    F2Counter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    F2Counter transversions;
    F2Counter transitions;
    //This excludes the following cases:
    // S = T, A = C
    // S = A, A = G
    F2Counter noDamage;

    F2Result();
    F2Result(const F2Result & other);
    ~F2Result();
    /* F2Result & operator= (const F2Result & other); */

    F2Result & operator= (const F2Result & other){
	this->all             = other.all;
	this->noCpg           = other.noCpg;
	this->onlyCpg         = other.onlyCpg;
	this->transversions   = other.transversions;
	this->transitions     = other.transitions;
	this->noDamage        = other.noDamage;
	return *this;
    }

    void addIfNotInfinity( vector<double> * target , double val );
    string printWithBootstrap(list<vector< vector< vector<F2Result> >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int k,unsigned int numberOfBootstraps);
    string printWithJacknife(const vector< const F2Result * >  * jacknife);
    void addAlleleCounts(int c1,int c2);

    F2Result &  operator+=(const F2Result & other);
    F2Result &  operator-=(const F2Result & other);

    friend istream &operator>>(istream &is , F2Result &dr){
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

    friend ostream& operator<<(ostream& os, const F2Result & dr){
	/* cout<<"F2Result print() begin"<<endl; */
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
