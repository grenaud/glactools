/*
 * AXT2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AXT2ACF_h
#define AXT2ACF_h

#include <string>
#include <climits>
#include <gzstream.h>

#include "libgab.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "GlactoolsOperations.h"

using namespace std;

class AXT2ACF{
private:
    bool uncompressed=0;    

public:
    AXT2ACF();
    AXT2ACF(const AXT2ACF & other);
    ~AXT2ACF();
    AXT2ACF & operator= (const AXT2ACF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
