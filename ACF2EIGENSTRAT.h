/*
 * Acf2eigenstrat
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2EIGENSTRAT_h
#define ACF2EIGENSTRAT_h

#include <string>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "utils.h"
#include "ReadTabix.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "GlactoolsOperations.h"

using namespace std;

class ACF2EIGENSTRAT{
private:

    bool printRoot=false;
    bool printAnc =false;
    bool limitToTransversions =false;

public:
    ACF2EIGENSTRAT();
    ACF2EIGENSTRAT(const ACF2EIGENSTRAT & other);
    ~ACF2EIGENSTRAT();
    ACF2EIGENSTRAT & operator= (const ACF2EIGENSTRAT & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
