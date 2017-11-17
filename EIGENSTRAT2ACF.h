/*
 * EIGENSTRAT2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef EIGENSTRAT2ACF_h
#define EIGENSTRAT2ACF_h

#include <string>
#include <climits>
#include <gzstream.h>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "GlactoolsOperations.h"

using namespace std;

class EIGENSTRAT2ACF{
private:
    bool uncompressed=0;    
    string     epoFile      = "none";
    bool       epoFileB     = false;
    string     fastaIndex   = "";
    void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO);

    const kstring_t * kstringPtrEPO;
    ks_tokaux_t aux;

public:
    EIGENSTRAT2ACF();
    EIGENSTRAT2ACF(const EIGENSTRAT2ACF & other);
    ~EIGENSTRAT2ACF();
    EIGENSTRAT2ACF & operator= (const EIGENSTRAT2ACF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
