/*
 * GlactoolsOperations
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlactoolsOperations_h
#define GlactoolsOperations_h

#include <vector>
#include "GlacParser.h"
/* #include <mlpack/core.hpp> */
#include "GlacWriter.h"
//#include "GlactoolsParser.h"
#include "GenomicRange.h"

/* using namespace arma; */
using namespace std;

typedef struct{
    string name;
    uint64_t startIndexChr;
    uint64_t endIndexChr;
    uint64_t length;
} chrinfo;

//reading fasta index
void readFastaIndex(const string fastaIndex,
		    vector<chrinfo> & chrFound,
		    uint64_t & genomeLength);

string initFiles(vector<GlacParser * > & vectorOfGP,
	       //bool & atLeastOneHasData,
		 vector<bool> & hasData,
		 vector<int> & popSizePerFile,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 //string & chr1,
		 uint16_t & chr1,
		 unsigned int & coordCurrent,
		 bool printOnlyFirstPop=false); //if we print only the populations of the first file

bool sanityCheck(vector<GlacParser * > & vectorOfGP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 //string & chr1,
		 uint16_t & chr1,
		 unsigned int & coordCurrent,
		 uint16_t & chrcheck,
		 //string & chrcheck ,
		 char & refAllele,
		 bool force=false);

bool printAllele(vector<GlacParser * > & vectorOfGP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<int> & popSizePerFile,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 //string & chr1,
		 uint16_t & chr1,
		 unsigned int & coordCurrent,
		 GlacWriter * gw,
		 bool isGL,
		 bool force=false);


map< string, vector<GenomicRange> * > * readBEDSortedfile(string filetoread,bool verbose=false);

inline int compare2ChrsU(const uint32_t & chr1,const uint32_t & chr2){
    // -1 if chr1<chr2
    if(chr1<chr2) return -1;
    //  0 if both chromosomes are equal
    if(chr1==chr2) return 0;
    //  1 if chr1>chr2
    if(chr1>chr2) return 1;

    return 0;//should not be reached
}




#endif
