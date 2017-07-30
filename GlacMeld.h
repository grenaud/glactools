/*
 * GlacCAT
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacCAT_h
#define GlacCAT_h

#include <string>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacCAT{
private:
    bool uncompressed=0;    
    int minPLdiffind=33;

public:
    GlacCAT();
    GlacCAT(const GlacCAT & other);
    ~GlacCAT();
    GlacCAT & operator= (const GlacCAT & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
