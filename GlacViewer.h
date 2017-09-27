/*
 * GlacViewer
 * Date: Jul-25-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacViewer_h
#define GlacViewer_h

#include <string>
#include <cinttypes>

#include "GlacWriter.h"


using namespace std;

class GlacViewer{
 private:
    bool printheader    =false;
    bool printdefline   =false;
    bool printBin       =false;
    bool uncompressed   =false;
    bool printpops      =false;
    bool printonlyheader=false;

    int bytesForAC=2;
    int bytesForGL=1;
    double subsampleProp=1.0;
    double subsampleB   =false;
    double headB        =false;
    unsigned int headBN =0;

 public:
    GlacViewer();
    GlacViewer(const GlacViewer & other);
    ~GlacViewer();
    GlacViewer & operator= (const GlacViewer & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
};


#endif
