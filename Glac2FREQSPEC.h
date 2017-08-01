/*
 * Glac2FREQSPEC
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Glac2FREQSPEC_h
#define Glac2FREQSPEC_h

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

class Glac2FREQSPEC{
private:

    bool onlysegsite = false;
    bool usefreq     = false;
    bool splitpop    = false;
    bool useAnc      = false;
    bool useRoot     = false;


public:
    Glac2FREQSPEC();
    Glac2FREQSPEC(const Glac2FREQSPEC & other);
    ~Glac2FREQSPEC();
    Glac2FREQSPEC & operator= (const Glac2FREQSPEC & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
