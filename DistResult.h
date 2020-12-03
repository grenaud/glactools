/*
 * DistResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef DistResult_h
#define DistResult_h

using namespace std;

#include <sstream>
#include <vector>
#include <list>
#include "DistAlleleCounter.h"
#include "libgab.h"

class DistResult{
 private:
 public:
    //This is everything
    DistAlleleCounter all;
    //This is without the ones marked as CpG
    DistAlleleCounter noCpg;
    //This is only with the ones marked as CpG
    DistAlleleCounter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    DistAlleleCounter transversions;


    DistAlleleCounter transitions;


    //This excludes the following cases:
    // S = T, R or A = C
    // S = A, R or A = G
    DistAlleleCounter noDamage;

    DistResult();
    ~DistResult();
    string getHeader();
    void addIfNotInfinity( vector<double> * target , double val );
    string printWithBootstrap(list< vector< vector< DistResult >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int numberOfBootstraps);
    string printWithJacknife(const vector< const DistResult * >  * jacknife,const string & dnaDistMode);
    string printAll() const;

    DistResult  & operator+= (const DistResult & other){
	all           += other.all;
	noCpg         += other.noCpg;
	onlyCpg       += other.onlyCpg;
	transitions   += other.transitions;
	transversions += other.transversions;
	
	return *this;
    }


    DistResult  & operator-= (const DistResult & other){
	all           -= other.all;
	noCpg         -= other.noCpg;
	onlyCpg       -= other.onlyCpg;
	transitions   -= other.transitions;
	transversions -= other.transversions;	
	return *this;
    }

    friend ostream& operator<<(ostream& os, const DistResult & ct){
	os<<ct.all<<"\t"
	  <<ct.onlyCpg<<"\t"
	  <<ct.noCpg<<"\t"
	  <<ct.transitions<<"\t"
	  <<ct.transversions<<"\t"
	  <<ct.noDamage;
	return os;
    }

    friend istream &operator>>(istream &is , DistResult &dr){
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

};
#endif
