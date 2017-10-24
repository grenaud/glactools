/*
 * FstResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef FstResult_h
#define FstResult_h

using namespace std;

#include <sstream>
#include <vector>
#include <list>
#include "FstAlleleCounter.h"
#include "utils.h"

class FstResult{
 private:
 public:
    //This is everything
    FstAlleleCounter all;
    //This is without the ones marked as CpG
    FstAlleleCounter noCpg;
    //This is only with the ones marked as CpG
    FstAlleleCounter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    FstAlleleCounter transversions;


    FstAlleleCounter transitions;


    //This excludes the following cases:
    // S = T, R or A = C
    // S = A, R or A = G
    FstAlleleCounter noDamage;

    FstResult();
    ~FstResult();
    string getHeader();
    void addIfNotInfinity( vector<double> * target , double val );
    string printWithBootstrap(list< vector< vector< FstResult >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int numberOfBootstraps);
    string printWithJacknife(const vector< const FstResult * >  * jacknife);

    FstResult  & operator+= (const FstResult & other){
	all           += other.all;
	noCpg         += other.noCpg;
	onlyCpg       += other.onlyCpg;
	transitions   += other.transitions;
	transversions += other.transversions;
	
	return *this;
    }


    FstResult  & operator-= (const FstResult & other){
	all           -= other.all;
	noCpg         -= other.noCpg;
	onlyCpg       -= other.onlyCpg;
	transitions   -= other.transitions;
	transversions -= other.transversions;	
	return *this;
    }

    friend ostream& operator<<(ostream& os, const FstResult & ct){
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
