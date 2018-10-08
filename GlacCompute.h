/*
 * GlacCompute
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacCompute_h
#define GlacCompute_h

#include <string>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <queue>
#include <map>
#include <algorithm>
#include <unistd.h>


#include "utils.h"
#include "ReadTabix.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "GlactoolsOperations.h"

#include "SumStatAvgCoa.h"
#include "SumStatFst.h"
#include "SumStatD.h"
#include "SumStatDist.h"
#include "SumStatF3.h"

using namespace std;

//! Chunk of code to check if a certain thread call failed
/*!
  This block is calls by the pthread

*/     
#define checkResults(string, val) {             \
 if (val) {					\
     cerr<<"Failed with "<<val<<" at "<<string<<endl;   \
   exit(1);                                     \
 }                                              \
}

class GlacCompute{
private:
    int numberOfThreads =       1;
    int sizeBins        = 1000000;
    bool performBoot=true;
    string program;
    string dnaDistMode;
public:
    GlacCompute();
    GlacCompute(const GlacCompute & other);
    ~GlacCompute();
    GlacCompute & operator= (const GlacCompute & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
