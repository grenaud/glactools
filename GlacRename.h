/*
 * GlacRename
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacRename_h
#define GlacRename_h

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

class GlacRename{
private:
    bool uncompressed=false;

public:
    GlacRename();
    GlacRename(const GlacRename & other);
    ~GlacRename();
    GlacRename & operator= (const GlacRename & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
