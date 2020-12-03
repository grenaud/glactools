/*
 * GlacPopsub
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacPopsub_h
#define GlacPopsub_h

#include <string>
#include <set>

#include "libgab.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacPopsub{
private:
    bool uncompressed=0;    

public:
    GlacPopsub();
    GlacPopsub(const GlacPopsub & other);
    ~GlacPopsub();
    GlacPopsub & operator= (const GlacPopsub & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
