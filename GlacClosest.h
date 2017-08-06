/*
 * GlacClosest
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacClosest_h
#define GlacClosest_h

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

class GlacClosest{
private:


public:
    GlacClosest();
    GlacClosest(const GlacClosest & other);
    ~GlacClosest();
    GlacClosest & operator= (const GlacClosest & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
