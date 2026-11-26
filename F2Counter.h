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

#include "libgab.h"

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
	double f2 = ct.computeF2();
	os<<ct.f2Sum<<"\t"<<ct.counterSites<<"\t"<<ct.nInds1<<"\t"<<ct.nInds2<<"\t"<<ct.hz1<<"\t"<<ct.hz2<<"\t"<<ct.nsnp<<"\t"<<f2<<"\t";
	return os;
    }

    friend istream &operator>>(istream &is , F2Counter &ct){
	is>>ct.f2Sum>>ct.counterSites>>ct.nInds1>>ct.nInds2>>ct.hz1>>ct.hz2>>ct.nsnp;

	string derived; //f2, recomputed from the sums
	is>>derived;

	return is;
    }


};



#endif
