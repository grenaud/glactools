/*
 * GlacUnion
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacUnion_h
#define GlacUnion_h

#include <string>

#include "libgab.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "GlactoolsOperations.h"

#include "AlleleRecords.h"

using namespace std;

class GlacUnion{
private:
    bool uncompressed=0;    
    bool force;
public:
    GlacUnion();
    GlacUnion(const GlacUnion & other);
    ~GlacUnion();
    GlacUnion & operator= (const GlacUnion & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
