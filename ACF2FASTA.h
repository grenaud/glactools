/*
 * ACF2FASTA
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ACF2FASTA_h
#define ACF2FASTA_h

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

class ACF2FASTA{
private:
    bool printRoot=true;
    bool produceTwo=false;


public:
    ACF2FASTA();
    ACF2FASTA(const ACF2FASTA & other);
    ~ACF2FASTA();
    ACF2FASTA & operator= (const ACF2FASTA & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
