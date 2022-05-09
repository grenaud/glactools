/*
 * BEAGLE2GLF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef BEAGLE2GLF_h
#define BEAGLE2GLF_h

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

#define MIN2(a,b) (((a)<(b))?(a):(b))
#define MIN3(a,b,c) MIN2(MIN2(a,b),c)

using namespace std;

class BEAGLE2GLF{
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
    BEAGLE2GLF();
    BEAGLE2GLF(const BEAGLE2GLF & other);
    ~BEAGLE2GLF();
    BEAGLE2GLF & operator= (const BEAGLE2GLF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
    inline int phredcapped(const double p) const;
};
#endif
