/*
 * RandomGenomicCoord
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef RandomGenomicCoord_h
#define RandomGenomicCoord_h

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <stdint.h>

#include "libgab.h"
#include "GenomicRange.h"
#include "GlactoolsOperations.h"


using namespace std;

/* typedef struct{ */
/*     string chrName; */
/*     unsigned int chrCoord; */
/* } randomCoord; */



/* typedef struct{ */
/*     string name; */
/*     uint64_t startIndexChr; */
/*     uint64_t endIndexChr; */
/*     uint64_t length; */
/* } chrinfoRandom; */


class RandomGenomicCoord{
 private:
    uint64_t genomeLength;
    vector<chrinfo> chrFound;
    uint64_t randGenomicCoord_();
    bool allowSexChr;
public:
    RandomGenomicCoord(string sqLines,bool allowSexChr=false );
    ~RandomGenomicCoord();

    GenomicRange getRandomGenomic(int bpToExtract);
};
#endif
