/*
 * GlacNosingle
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacNosingle_h
#define GlacNosingle_h

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

class GlacNosingle{
private:
    bool uncompressed  =false;
    bool includeRootAnc=false;

public:
    GlacNosingle();
    GlacNosingle(const GlacNosingle & other);
    ~GlacNosingle();
    GlacNosingle & operator= (const GlacNosingle & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
