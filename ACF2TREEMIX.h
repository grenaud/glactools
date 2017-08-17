/*
 * ACF2TREEMIX
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2TREEMIX_h
#define ACF2TREEMIX_h

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

class ACF2TREEMIX{
private:

    bool limitToTransversions=false;
    bool noprivate=false;
    bool printAnc=false;

public:
    ACF2TREEMIX();
    ACF2TREEMIX(const ACF2TREEMIX & other);
    ~ACF2TREEMIX();
    ACF2TREEMIX & operator= (const ACF2TREEMIX & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
