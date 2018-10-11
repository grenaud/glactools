/*
 * GlacWindows
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacWindows_h
#define GlacWindows_h

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
#include "GenomicWindows.h"
#include "GlactoolsOperations.h"

using namespace std;

class GlacWindows{
private:
    //    bool inBedFormat=false;
    // case 'v':  //entire chr (only one window)
    char lociSelection='x';

    int bpToExtract =-1      ; //amount of bp per chunk
    int overlapBetweenWindows=0;
    string chrToUse;
    unsigned int amountOfGoodSitesTARGET=0;  //number of contiguous chunks to extract for random
    bool allowSexChr=false;

    bool specifiedChunk=false;
    bool specifiedSizeFile=false;
    string sizeFile;
    bool inBedFormat=false;

public:
    GlacWindows();
    GlacWindows(const GlacWindows & other);
    ~GlacWindows();
    GlacWindows & operator= (const GlacWindows & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
