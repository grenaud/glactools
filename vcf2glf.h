/*
 * vcf2glf
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef vcf2glf_h
#define vcf2glf_h


#include "utils.h"
#include "SetVCFFilters.h"
#include "SimpleVCF.h"
#include "VCFreader.h"
#include "AlleleInfo.h"
#include "FilterVCF.h"
#include "glactoolsOperations.h"
#include <cinttypes>

extern "C"{
#include "htslib/bgzf.h"
}


using namespace std;

class vcf2glf{
private:

int minPLdiffind=33;
int minGQcutoff=0;
int minMQcutoff=0;
int bytesForGL=1;
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
public:
    vcf2glf();
    vcf2glf(const vcf2glf & other);
    ~vcf2glf();
    vcf2glf & operator= (const vcf2glf & other);
    
    string usage() const;
    int run(int argc, char *argv[]);
};
#endif
