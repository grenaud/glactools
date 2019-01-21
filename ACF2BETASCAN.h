/*
 * ACF2BETASCAN
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2BETASCAN_h
#define ACF2BETASCAN_h

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

class ACF2BETASCAN{
private:

    /* bool onlysegsite = false; */
    /* bool usefreq     = false; */
    /* bool splitpop    = false; */

    bool fold        = false;
    bool useAnc      = false;
    bool useRoot     = false;


public:
    ACF2BETASCAN();
    ACF2BETASCAN(const ACF2BETASCAN & other);
    ~ACF2BETASCAN();
    ACF2BETASCAN & operator= (const ACF2BETASCAN & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
