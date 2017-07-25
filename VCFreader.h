/*
 * VCFreader
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef VCFreader_h
#define VCFreader_h

#include <list>
#include <fstream> 
#include <memory> 
#include <gzstream.h>

#include "SimpleVCF.h"
#include "ReadTabix.h"
#include "AlleleInfoReader.h"

using namespace std;

//! This class is to read a simple VCF (not multi individuals)
/*!
 *
 *  This class will read VCF data and flag those within a certain 
 *  range (ex: 5bp) as close to an indel.
 *
 *
 *    VCFreader vcfr ("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz",
 *		    "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz.tbi",
 *		    "10",
 *		    557572,
 *		    557592,
 *		    5);
 *
 *		    // 75060,
 *		    // 75070,
 *		    // 5);
 *    
 *    while(vcfr.hasData()){
 *	SimpleVCF toprint=vcfr.getData();
 *	cout<<toprint<<"\t"<<toprint.containsIndel()<<"\t"<<toprint.getCloseIndel()<<endl;
 *    }
 *
 *
 *
 *    VCFreader vcfr2 (string(argv[1]),
 *		     7);
 *
 *    
 *    while(vcfr2.hasData()){
 *	SimpleVCF toprint=vcfr2.getData();
 *	cout<<toprint<<"\t"<<toprint.containsIndel()<<"\t"<<toprint.getCloseIndel()<<endl;
 *    }
 *
 *
*/
class VCFreader : public AlleleInfoReader{
private:
    //childproof variable to make sure the user called hasData() before calling getData()
    int numberOfTimesHasDataWasCalled;

    ReadTabix * rt;
    string currentline;
    int readAhead;
    list<SimpleVCF *> queueOfVCFs;
    int indexInQueueOfIndels; //this is the index of the last element of the list that is an insert, -1 if none

    //data of last indel found
    string chrOfLastIndel;
    unsigned int indexOfLastIndel;
    bool previouslyFoundIndel;

    bool needToPopulateQueue;
    bool fullQueue;
    bool endQueue;

    bool tabixMode;
    bool textMode;

    bool repoCalledHasData; //flag if the repositionIterator was called, set to true by repositionIterator() and false by hasData()

    igzstream vcfFile; //for text mode
    inline bool getNextLine();
    SimpleVCF * svcfToReturn;
    inline void flagCpG(SimpleVCF * previous,SimpleVCF * current);

public:


//! To retrieve the records with tabix
/*!
 *
 * This constructor retrieves a section of a file using tabix

  \param file : Full path of the vcf file (zipped with bgzip beforehand)
  \param indexForFile : Full path of the tabix index for the  vcf file
  \param chrName : Chromosome name
  \param start : Start coordinate 
  \param end : end coordinate 
  \param indelsAhead : Mark this many positions around an indel to flag them as close to indel
*/
    VCFreader(string file,string indexForFile,string chrName,int start,int end,int indelsAhead=5);


//! To retrieve the records for a simple uncompressed VCF file
/*!
 *
 * This constructor retrieves a every VCF record in an uncompressed VCF file

  \param file : Full path of the vcf file
  \param indelsAhead : Mark this many positions around an indel to flag them as close to indel
*/
    VCFreader(string file,int indelsAhead=5);



    ~VCFreader();



//! To reposition the tabix iterator
/*!
 *
 * This subroutine repositions the tabix iterator if it was opened using tabix (on the same file)
 * for a diverent file or change the number of indels, call the constructor again

  \param chrName : Chromosome name
  \param start : Start coordinate 
  \param end : end coordinate 
*/
    void repositionIterator(string chrName,int start,int end);



//! To verify if there are still data to be read
/*!
 *
 * This subroutine must be called PRIOR to calling getData
 *
  \return true if there is still data
  \sa getData
*/
    bool hasData();

//! To retrieve the data from the next line
/*!
 *
 * This subroutine must be called AFTER calling hasData
 *
  \return a SimpleVCF object
  \sa hasData
*/
    SimpleVCF * getData();
    /* auto_ptr<AlleleInfo> getData(); */

};
#endif
