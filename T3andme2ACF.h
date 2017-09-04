/*
 * T3andme2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef T3andme2ACF_h
#define T3andme2ACF_h

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

class T3andme2ACF{
private:
bool uncompressed=0;    
string epoFile   = "none";
bool   epoFileB  = false;
string     fastaIndex   =""  ;

public:
    T3andme2ACF();
    T3andme2ACF(const T3andme2ACF & other);
    ~T3andme2ACF();
    T3andme2ACF & operator= (const T3andme2ACF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
    void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO);
};
#endif
