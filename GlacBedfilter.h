/*
 * GlacBedfilter
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacBedfilter_h
#define GlacBedfilter_h

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

class GlacBedfilter{
private:
    bool verbose=false;
    bool uncompressed=false;

public:
    GlacBedfilter();
    GlacBedfilter(const GlacBedfilter & other);
    ~GlacBedfilter();
    GlacBedfilter & operator= (const GlacBedfilter & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
