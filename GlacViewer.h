/*
 * GlacViewer
 * Date: Jul-25-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacViewer_h
#define GlacViewer_h
#include <string>

using namespace std;

class GlacViewer{
 private:
    bool printheader =false;
    bool printdefline=false;
    bool printBin      =false;
    bool uncompressed  =false;
    int bytesForAC=2;
    int bytesForGL=1;

 public:
    GlacViewer();
    GlacViewer(const GlacViewer & other);
    ~GlacViewer();
    GlacViewer & operator= (const GlacViewer & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
};


#endif
