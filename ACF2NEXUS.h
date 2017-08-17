/*
 * ACF2NEXUS
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2NEXUS_h
#define ACF2NEXUS_h

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

class ACF2NEXUS{
private:
    bool printOnlySeq =false;
    bool allDefined=false;


public:
    ACF2NEXUS();
    ACF2NEXUS(const ACF2NEXUS & other);
    ~ACF2NEXUS();
    ACF2NEXUS & operator= (const ACF2NEXUS & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
