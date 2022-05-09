/*
 * GlacUndef
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacUndef_h
#define GlacUndef_h

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

class GlacUndef{
private:
    bool uncompressed=false;

public:
    GlacUndef();
    GlacUndef(const GlacUndef & other);
    ~GlacUndef();
    GlacUndef & operator= (const GlacUndef & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
