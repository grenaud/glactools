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

#include "libgab.h"
#include "CheckPrintedValue.h"

//#include <iostream>
//#include <fstream>
using namespace std;

class F3Counter{
 private:

 public:

    //counter of mutations

    double   f3Sum;
    unsigned int counterSites;
    double nInds; 
    double hzy; 
    double nsnp; 


    F3Counter(); //constructor
    //    vector<double>   performBoostraps() const;
    void reinitializedCounters();
    string headerForCount() const;

    double computeF3() const;
    void addAlleleCounts(int c1,int c2);

    F3Counter &  operator+=(const F3Counter & other);
    F3Counter &  operator-=(const F3Counter & other);

    friend ostream& operator<<(ostream& os, const F3Counter & ct){
	double f3unscaled = double(ct.f3Sum)/double(ct.counterSites);
	double t          = double(ct.hzy)/double(ct.counterSites);//mean correction term

	os<<ct.f3Sum<<"\t"<<ct.counterSites<<"\t"<<f3unscaled<<"\t"<<ct.nInds<<"\t"<<ct.hzy<<"\t"<<ct.nsnp<<"\t"<<t<<"\t"<<ct.computeF3()<<"\t";
	return os;
    }

    friend istream &operator>>(istream &is , F3Counter &ct){
	//the derived columns are recomputed from the sums
	string derived,printed;
	is>>ct.f3Sum>>ct.counterSites>>derived;      //uncorrected f3
	is>>ct.nInds>>ct.hzy>>ct.nsnp>>derived>>printed; //correction term, corrected f3
	checkPrintedValue("f3",printed,ct.computeF3(),(std::fabs(ct.f3Sum)+ct.hzy)/double(ct.counterSites));

	return is;
    }


};



#endif
