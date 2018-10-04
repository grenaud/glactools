/*
 * F3Counter
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef F3Counter_h
#define F3Counter_h

#include <iostream>
#include <ostream> 
#include <vector> 
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

#include "utils.h"

//#include <iostream>
//#include <fstream>
using namespace std;

class F3Counter{
 private:

 public:

    //counter of mutations

    double   f3Sum;
    unsigned int counterSites;


    F3Counter(); //constructor
    vector<double>   performBoostraps() const;
    void reinitializedCounters();
    string headerForCount() const;

    double computeDST() const;

    F3Counter &  operator+=(const F3Counter & other);
    F3Counter &  operator-=(const F3Counter & other);

    friend ostream& operator<<(ostream& os, const F3Counter & ct){
	/* vector<double> boostraps = ct.performBoostraps(); */
	//TODO
	/* double dst= (double(ct.counterAncDer)-double(ct.counterDerAnc))/(double(ct.counterAncDer)+double(ct.counterDerAnc));  */
	/* os<<ct.counterAncAnc<<"\t" */
	/*   <<ct.counterAncDer<<"\t" */
	/*   <<ct.counterDerAnc<<"\t" */
	/*   <<ct.counterDerDer<<"\t" */
	/*     //(ADDA - DADA )/ (ADDA+DADA ) */
	/*   <<dst<<"\t" ; */
	/* pair<double,double> res = computeMeanSTDDEV(boostraps); */
	/* os<< ((res.first-dst)/res.second); */
	return os;
    }


};



#endif
