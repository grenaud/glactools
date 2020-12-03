/*
 * GlacNoStrictSharing
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacNoStrictSharing_h
#define GlacNoStrictSharing_h

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

class GlacNoStrictSharing{
private:
    bool uncompressed=false;

public:
    GlacNoStrictSharing();
    GlacNoStrictSharing(const GlacNoStrictSharing & other);
    ~GlacNoStrictSharing();
    GlacNoStrictSharing & operator= (const GlacNoStrictSharing & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
