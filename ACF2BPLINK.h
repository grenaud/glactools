/*
 * Acf2bplink
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2BPLINK_h
#define ACF2BPLINK_h

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

class ACF2BPLINK{
private:

    bool printRoot    = true;
    bool singleAlHomo = false;
    bool haploidRoot  = false;

public:
    ACF2BPLINK();
    ACF2BPLINK(const ACF2BPLINK & other);
    ~ACF2BPLINK();
    ACF2BPLINK & operator= (const ACF2BPLINK & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
