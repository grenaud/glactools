/*
 * GlacDown
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacDown_h
#define GlacDown_h

#include <string>
#include <set>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacDown{
private:
    bool uncompressed=0;    
    int count        =1;

public:
    GlacDown();
    GlacDown(const GlacDown & other);
    ~GlacDown();
    GlacDown & operator= (const GlacDown & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
