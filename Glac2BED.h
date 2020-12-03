/*
 * Glac2BED
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Glac2BED_h
#define Glac2BED_h

#include <string>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "libgab.h"
#include "ReadTabix.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "GlactoolsOperations.h"

using namespace std;

class Glac2BED{
private:


public:
    Glac2BED();
    Glac2BED(const Glac2BED & other);
    ~Glac2BED();
    Glac2BED & operator= (const Glac2BED & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
