/*
 * Glac2VCF
 * Date: Apr-10-2019
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GLAC2VCF_h
#define GLAC2VCF_h

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

class Glac2VCF{
private:

    bool printRoot          = false;
    bool printAnc           = false;
    bool singleAlleleAsHomo = false;

public:
    Glac2VCF();
    Glac2VCF(const Glac2VCF & other);
    ~Glac2VCF();
    Glac2VCF & operator= (const Glac2VCF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
