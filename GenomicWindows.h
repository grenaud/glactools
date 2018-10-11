/*
 * GenomicWindows
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef GenomicWindows_h
#define GenomicWindows_h

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "utils.h"
#include "GenomicRange.h"
#include "RandomGenomicCoord.h"


using namespace std;



class GenomicWindows{
 private:
    unsigned int genomeLength;
    vector<chrinfo> chrFound;
    bool allowSexChr;

public:
    GenomicWindows(string sqLines,bool allowSexChr=false );
    ~GenomicWindows();

    vector<GenomicRange> getGenomicWindows(int windowSize,int overlap=0);
    vector<GenomicRange> getGenomicWindowsChr(string chrName, int windowSize,int overlap=0);


    vector<GenomicRange> getGenomeWide();
    vector<GenomicRange> getChr(string chrName);




};
#endif
