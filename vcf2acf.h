/*
 * vcf2acf
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef vcf2acf_h
#define vcf2acf_h


#include "utils.h"
#include "SetVCFFilters.h"


using namespace std;

class vcf2acf{
private:
    int minPLdiffind=33;
    int minGQcutoff=0;
    int minMQcutoff=0;

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

public:
    vcf2acf();
    vcf2acf(const vcf2acf & other);
    ~vcf2acf();
    vcf2acf & operator= (const vcf2acf & other);
    
    string usage() const;
    //void run(int argc, char *argv[]);
};
#endif
