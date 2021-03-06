/*
 * VcfGLF
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Vcf2GLF_h
#define Vcf2GLF_h


#include "libgab.h"
#include "SetVCFFilters.h"
#include "SimpleVCF.h"
#include "VCFreader.h"
#include "AlleleInfo.h"
#include "FilterVCF.h"
#include "GlacWriter.h"
#include "GlactoolsOperations.h"
#include <cinttypes>
#include <climits>

/* extern "C"{ */
/* #include "htslib/bgzf.h" */
/* } */


using namespace std;

class Vcf2GLF{
private:
    string epoFile   = "none";
    bool   epoFileB  = false;

int minPLdiffind=33;
int minGQcutoff=0;
int minMQcutoff=0;
int bytesForAC=2;
int minCovcutoff=0;
int maxCovcutoff=1000;
SetVCFFilters  * filtersVCF;
bool allowCloseIndelProx = false;
bool allowRepeatMasking  = false;
bool allowSysErr         = false;
bool allowall            = false;
bool allowallMQ          = false;

// bool ancAllele           = false;

double     minMapabilitycutoff =0;
string     fastaIndex   =""  ;
bool uncompressed=0;    
    void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO);

public:
    Vcf2GLF();
    Vcf2GLF(const Vcf2GLF & other);
    ~Vcf2GLF();
    Vcf2GLF & operator= (const Vcf2GLF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
};
#endif
