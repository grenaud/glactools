/*
 * GlacIndex
 * Date: Jul-25-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacIndex_h
#define GlacIndex_h

#include <iostream>
#include <vector>

#include <string>

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#include "hts_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "utils.h"

using namespace std;


class GlacIndex{
private:

public:
    GlacIndex();
    GlacIndex(const GlacIndex & other);
    ~GlacIndex();
    GlacIndex & operator= (const GlacIndex & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
    
};
#endif
