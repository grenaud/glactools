/*
 * GlacUsePopAsRootAnc
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacUsePopAsRootAnc_h
#define GlacUsePopAsRootAnc_h

#include <string>
#include <set>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacUsePopAsRootAnc{
private:
    bool uncompressed=0;    

public:
    GlacUsePopAsRootAnc();
    GlacUsePopAsRootAnc(const GlacUsePopAsRootAnc & other);
    ~GlacUsePopAsRootAnc();
    GlacUsePopAsRootAnc & operator= (const GlacUsePopAsRootAnc & other);
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
