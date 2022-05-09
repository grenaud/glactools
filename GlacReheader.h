/*
 * GlacReheader
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacReheader_h
#define GlacReheader_h

#include <string>
#include <set>

#include "libgab.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacReheader{
private:
    bool uncompressed=0;    

public:
    GlacReheader();
    GlacReheader(const GlacReheader & other);
    ~GlacReheader();
    GlacReheader & operator= (const GlacReheader & other);
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
