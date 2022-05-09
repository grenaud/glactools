/*
 * GlacStats
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacStats_h
#define GlacStats_h

#include <string>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "libgab.h"
#include "ReadTabix.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "GlactoolsOperations.h"

using namespace std;

class GlacStats{
private:

    /* bool onlysegsite = false; */
    /* bool usefreq     = false; */
    /* bool splitpop    = false; */
    /* bool useAnc      = false; */
    /* bool useRoot     = false; */


public:
    GlacStats();
    GlacStats(const GlacStats & other);
    ~GlacStats();
    GlacStats & operator= (const GlacStats & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
