/*
 * AvgCoaResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AvgCoaResult_h
#define AvgCoaResult_h

using namespace std;

#include <sstream>
#include <vector>
#include <list>
#include "AvgCoaAlleleCounter.h"
#include "utils.h"

class AvgCoaResult{
 private:
 public:
    //This is everything
    AvgCoaAlleleCounter all;
    //This is without the ones marked as CpG
    AvgCoaAlleleCounter noCpg;
    //This is only with the ones marked as CpG
    AvgCoaAlleleCounter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    AvgCoaAlleleCounter transversions;


    AvgCoaAlleleCounter transitions;


    //This excludes the following cases:
    // S = T, R or A = C
    // S = A, R or A = G
    AvgCoaAlleleCounter noDamage;

    AvgCoaResult();
    ~AvgCoaResult();
    string getHeader();
    void addIfNotInfinity( vector<double> * target , double val );
    string printWithBootstrap(list< vector< vector< AvgCoaResult >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int numberOfBootstraps);
    string printWithJacknife(const vector< const AvgCoaResult * >  * jacknife);

    AvgCoaResult  & operator+= (const AvgCoaResult & other){
	all           += other.all;
	noCpg         += other.noCpg;
	onlyCpg       += other.onlyCpg;
	transitions   += other.transitions;
	transversions += other.transversions;
	
	return *this;
    }


    AvgCoaResult  & operator-= (const AvgCoaResult & other){
	all           -= other.all;
	noCpg         -= other.noCpg;
	onlyCpg       -= other.onlyCpg;
	transitions   -= other.transitions;
	transversions -= other.transversions;	
	return *this;
    }

    friend istream &operator>>(istream &is , AvgCoaResult &dr){
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

    friend ostream& operator<<(ostream& os, const AvgCoaResult & ct){
	os<<ct.all<<"\t"
	  <<ct.onlyCpg<<"\t"
	  <<ct.noCpg<<"\t"
	  <<ct.transitions<<"\t"
	  <<ct.transversions<<"\t"
	  <<ct.noDamage;
	return os;
    }
};
#endif
