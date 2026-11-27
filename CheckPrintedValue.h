/*
 * CheckPrintedValue
 *
 */

#ifndef CheckPrintedValue_h
#define CheckPrintedValue_h

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

//For --boot: the statistic printed on a --noboot row must agree with the sums on that row. It does not if
//the file was written by a version that stored other quantities in the same columns.
//magnitude is the size of the terms the statistic is the difference of.
inline void checkPrintedValue(const string & statistic,const string & printed,double recomputed,double magnitude){
    char * end=NULL;
    double value = strtod(printed.c_str(),&end);//also reads nan, -nan and inf
    if(!printed.empty() && *end == '\0'){
	if(std::isnan(value) && std::isnan(recomputed))
	    return;
	//both are computed from numbers printed with 6 significant digits
	if(std::fabs(value-recomputed) <= 2e-5*(magnitude+std::fabs(value)) + 1e-9)
	    return;
    }
    cerr<<"ERROR: the "<<statistic<<" value \""<<printed<<"\" in the input does not match the sums on the same row. The file was probably written by glactools 1.0.10 or earlier, which stored the sums differently: run \"glactools compute --noboot\" again with this version"<<endl;
    exit(1);
}

#endif
