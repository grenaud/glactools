/*
 * ACF2GPHOCS
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2GPHOCS_h
#define ACF2GPHOCS_h

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

class ACF2GPHOCS{
private:
    bool allowCpg=false;

public:
    ACF2GPHOCS();
    ACF2GPHOCS(const ACF2GPHOCS & other);
    ~ACF2GPHOCS();
    ACF2GPHOCS & operator= (const ACF2GPHOCS & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
