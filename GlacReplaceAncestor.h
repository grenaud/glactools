/*
 * GlacReplaceAncestor
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacReplaceAncestor_h
#define GlacReplaceAncestor_h

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

class GlacReplaceAncestor{
private:
    bool uncompressed=false;

public:
    GlacReplaceAncestor();
    GlacReplaceAncestor(const GlacReplaceAncestor & other);
    ~GlacReplaceAncestor();
    GlacReplaceAncestor & operator= (const GlacReplaceAncestor & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
