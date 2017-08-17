/*
 * GlacSharing
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacSharing_h
#define GlacSharing_h

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

class GlacSharing{
private:
    bool uncompressed=false;

public:
    GlacSharing();
    GlacSharing(const GlacSharing & other);
    ~GlacSharing();
    GlacSharing & operator= (const GlacSharing & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
