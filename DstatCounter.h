/*
 * DstatCounter
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef DstatCounter_h
#define DstatCounter_h

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

class DstatCounter{
 private:

 public:

    //counter of mutations
    unsigned int counterAncAnc;
    unsigned int counterAncDer;
    unsigned int counterDerAnc;
    unsigned int counterDerDer;


    DstatCounter(); //constructor
    //    vector<double>   performBoostraps() const;
    void reinitializedCounters();
    string headerForCount() const;

    double computeDST() const;

    DstatCounter &  operator+=(const DstatCounter & other);
    DstatCounter &  operator-=(const DstatCounter & other);

    friend ostream& operator<<(ostream& os, const DstatCounter & ct){
	/* vector<double> boostraps = ct.performBoostraps(); */
	double dst= (double(ct.counterAncDer)-double(ct.counterDerAnc))/(double(ct.counterAncDer)+double(ct.counterDerAnc)); 
	os<<ct.counterAncAnc<<"\t"
	  <<ct.counterAncDer<<"\t"
	  <<ct.counterDerAnc<<"\t"
	  <<ct.counterDerDer<<"\t"
	    //(ADDA - DADA )/ (ADDA+DADA )
	  <<dst<<"\t" ;
	/* pair<double,double> res = computeMeanSTDDEV(boostraps); */
	/* os<< ((res.first-dst)/res.second); */
	return os;
    }


};



#endif
