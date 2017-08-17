/*
 * GlacMeld
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacMeld_h
#define GlacMeld_h

#include <string>
#include <set>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacMeld{
private:
    bool uncompressed=0;    
    //int minPLdiffind=33;
    bool keepOrig=false;

public:
    GlacMeld();
    GlacMeld(const GlacMeld & other);
    ~GlacMeld();
    GlacMeld & operator= (const GlacMeld & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
