/*
 * BAM2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef BAM2ACF_h
#define BAM2ACF_h

#include <string>


#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>

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

};
#endif
