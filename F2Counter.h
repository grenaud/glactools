/*
 * F2Counter
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef F2Counter_h
#define F2Counter_h

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

class F2Counter{
 private:

 public:

    //counter of mutations

    double   f2Sum; //(f1-f2)^2
    unsigned int counterSites;

    double nInds1; 
    double nInds2; 

    double hz1;
    double hz2;
 
    /* double nsnp1;  */
    /* double nsnp2;  */
    double nsnp; 



    F2Counter(); //constructor
    //    vector<double>   performBoostraps() const;
    void reinitializedCounters();
    string headerForCount() const;

    double computeF2() const;
    void addAlleleCounts(int c1_1,int c2_1,double f_1,int c1_2,int c2_2,double f_2);

    F2Counter &  operator+=(const F2Counter & other);
    F2Counter &  operator-=(const F2Counter & other);


	/* double toreturn = 0; */
	/* for (int i = 0; i < nsnp; i++){ */
	/* 	double n1 = mean_ninds[p1]; */
	/* 	double n2 = mean_ninds[p2]; */
	/* 	double f1 = gsl_matrix_get(alfreqs, i, p1); */
	/* 	double f2 = gsl_matrix_get(alfreqs, i, p2); */

	/* 	double diff = f1-f2; */
	/* 	double hz1 = f1*(1-f1); */
	/* 	double hz2 = f2*(1-f2); */
	/* 	double toadd = diff*diff - hz1/(2*n1) - hz2/(2*n2); */
	/* 	toreturn += toadd; */
	/* } */
	/* toreturn = toreturn/ (double) nsnp; */
	/* return toreturn; */
    //taken from countData from treemix
    friend ostream& operator<<(ostream& os, const F2Counter & ct){
	double f2unscaled = double(ct.f2Sum);///double(ct.counterSites);

	double meanNind1   = double(ct.nInds1)/double(ct.nsnp);
	double meanNind2   = double(ct.nInds2)/double(ct.nsnp);

	double meanHZ1     = double(ct.hz1)/double(2*meanNind1);
	double meanHZ2     = double(ct.hz2)/double(2*meanNind2);

	double f2=  ( f2unscaled - meanHZ1 - meanHZ2 )  / double(ct.nsnp);
	//  324.75	  9344	    324.75	9304.5	9335	227.5	217.25	9344	0.0108933
	//  324.75	         9344	               9304.5	          9335	         227.5         217.25	     9344	0.0108933
	os<<ct.f2Sum<<"\t"<<ct.counterSites<<"\t"<<ct.nInds1<<"\t"<<ct.nInds2<<"\t"<<ct.hz1<<"\t"<<ct.hz2<<"\t"<<ct.nsnp<<"\t"<<f2<<"\t";
	return os;
    }


};



#endif
