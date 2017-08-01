/*
 * GlacIDXSTATS
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacIDXSTATS_h
#define GlacIDXSTATS_h

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

class GlacIDXSTATS{
private:


public:
    GlacIDXSTATS();
    GlacIDXSTATS(const GlacIDXSTATS & other);
    ~GlacIDXSTATS();
    GlacIDXSTATS & operator= (const GlacIDXSTATS & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
