/*
 * GlacSegsite
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacSegsite_h
#define GlacSegsite_h

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

class GlacSegsite{
private:
    bool uncompressed=false;

    bool onlyTS=false;
    bool onlyTV=false;

public:
    GlacSegsite();
    GlacSegsite(const GlacSegsite & other);
    ~GlacSegsite();
    GlacSegsite & operator= (const GlacSegsite & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
