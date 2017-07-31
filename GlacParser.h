#ifndef GlacParser_h
#define GlacParser_h

//#include <gzstream.h>
#include <sys/time.h> //for srand
#include <climits>

#include "utils.h"
#include "SingleAllele.h"
#include "AlleleRecords.h"
#include "ReadTabix.h"

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#include "hts_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/kstring.h"


using namespace std;


static uint8_t maskRef=240; //11110000
static uint8_t maskAlt=15; //00001111

static char sizeBytesACF;
static char sizeBytesGLF;


class GlacParser{
 private:

    vector<string> * populationNames;
    unsigned int numberPopulations;
    //igzstream   * myFilezipped;
    BGZF *myFilezipped; 
hts_itr_t *iter;//for iterator for indexing

    /* ifstream    * myFile; */
    AlleleRecords * allRecToReturn;
    int numberOfTimesHasDataWasCalled;
    string header;
    string headerNoSQNoDefline;
    string headerSQ;
    string headerNoDefline;
    const vector<string> * dataToRead;
    unsigned int dataToReadInd;

    string defline;
    string currentline;

    int numbernew;
    int numberdel;

    ReadTabix * rt;

    bool stringMode;
    bool tabixMode;
    bool textMode;
    bool binMode;
    bool acFormat;
    bool glFormat;

    uint32_t sizePops;

    vector<string> chrKnown;
    //void parseHeader(istream & in);
    void parseHeader(BGZF *myFilezipped); 
    /* bool getNextLine(); */

 public:
    GlacParser(string filename,int compressionThreads=1);
    GlacParser(string file,string indexForFile,string chrName,int start,int end,bool justChr=false,int compressionThreads=1);
    /* GlacParser(string file,string indexForFile); */

    //GlacParser(const vector<string> * dataToRead,const vector<string> & populationNames_);
   
    GlacParser(const GlacParser&); // not implemented
    ~GlacParser();
    bool hasData();
    AlleleRecords  * getData();

    string getHeader(string prefix="") const;
    string getHeaderNoSQNoDefline(string prefix="") const;
    string getHeaderSQ(string prefix="") const;
    string getHeaderNoDefline(string prefix="") const;

    string getDefline() const;
/* void repositionIterator(string chrName,int start,int end); */
    bool isACFormat() const;
    bool isGLFormat() const;
    uint32_t getSizePops() const;
    size_t getSizeRecord() const; //size of 1 record in binary

    const vector<string> *   getPopulationsNames() const ;

};




inline char sampleRandomRefAltAllele(char ref,char alt,int refCount,int altCount);

inline char sampleRandomRefAltAllele(char ref,char alt,int refCount,int altCount){

    if(refCount == 0 && altCount == 0 ){
	cerr<<"GlacParser::sampleRandomRefAltAllele() cannot sample when both allele counts are 0, find the null record in your input, exiting"<<endl;
	exit(1);
    }

    if(refCount == 0  ){//no reference, return alternative
	return alt;
    }

    if(altCount == 0  ){//no reference, return alternative
	return ref;
    }

    if(refCount ==  altCount ){//homozygous (kind of)
	if(randomBool())
	    return ref;
	else
	    return alt;
    }

    int totalCount = refCount+  altCount;
    int randnum=(callRand()%totalCount)+1;
    /* cout<<refCount<<" "<<altCount<<" randnum "<<randnum<<endl; */

    if(randnum<=refCount){
	/* cout<<"REF"<<endl; */
	return ref;
    }else{
	/* cout<<"ALT"<<endl; */
	return alt;
    }
    
}

#endif
