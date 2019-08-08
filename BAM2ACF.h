/*
 * BAM2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef BAM2ACF_h
#define BAM2ACF_h

#include <string>


extern "C" {
    //#include "tabix.h"
    //#include "bam.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "bam.h"

#include "samtools.h"
#include "sam_opts.h"
#include "bedidx.h"
}

#define bam_is_properpair(b)    (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam_is_paired(b)    (((b)->core.flag&BAM_FPAIRED) != 0)
#define bam_is_pair1(b)     (((b)->core.flag&BAM_FREAD1)  != 0)
#define bam_is_pair2(b)     (((b)->core.flag&BAM_FREAD2)  != 0)
#define bam_is_qcfailed(b)  (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_rmdup(b)     (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_failed(b)    ( bam_is_qcfailed(b) || bam_is_rmdup(b) )
#define bam_mqual(b)        ((b)->core.qual)
#define bam1_qname(b)       (bam_get_qname((b)))

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

/* #include <api/BamConstants.h> */
/* #include <api/BamMultiReader.h> */

#include <utils/bamtools_fasta.h>

/* #include <utils/bamtools_options.h> */
/* #include <utils/bamtools_pileup_engine.h> */
/* #include <utils/bamtools_utilities.h> */

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

class BAM2ACF{
private:

    bool uncompressed =  0;    

    string epoFile    = "none";
    bool   epoFileB   = false;
    bool useQCFail    = false;
    int minBaseQual   = 0;

public:
    BAM2ACF();
    BAM2ACF(const BAM2ACF & other);
    ~BAM2ACF();
    BAM2ACF & operator= (const BAM2ACF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
    static void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO);
};
#endif
