/*
 * GlacNoSharing
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacNoSharing_h
#define GlacNoSharing_h

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

class GlacNoSharing{
private:
    bool uncompressed=false;

public:
    GlacNoSharing();
    GlacNoSharing(const GlacNoSharing & other);
    ~GlacNoSharing();
    GlacNoSharing & operator= (const GlacNoSharing & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
