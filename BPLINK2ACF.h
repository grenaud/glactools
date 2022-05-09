/*
 * BPLINK2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef BPLINK2ACF_h
#define BPLINK2ACF_h

#include <string>
#include <climits>
#include <gzstream.h>

/* #include <api/BamConstants.h> */
/* #include <api/BamMultiReader.h> */
#include <utils/bamtools_fasta.h>
/* #include <utils/bamtools_options.h> */
/* #include <utils/bamtools_utilities.h> */


#include "libgab.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "GlactoolsOperations.h"

using namespace std;

class BPLINK2ACF{
private:
    bool uncompressed=0;    
    string     epoFile      = "none";
    bool       epoFileB     = false;
    string     fastaFile    = "";
    string     fastaIndex   = "";
    void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO);

    const kstring_t * kstringPtrEPO;
    ks_tokaux_t aux;

public:
    BPLINK2ACF();
    BPLINK2ACF(const BPLINK2ACF & other);
    ~BPLINK2ACF();
    BPLINK2ACF & operator= (const BPLINK2ACF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
