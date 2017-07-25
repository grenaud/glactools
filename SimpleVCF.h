/*
 * SimpleVCF
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef SimpleVCF_h
#define SimpleVCF_h

#include <string> 
#include <vector> 
#include <map>

#include "CoreVCF.h"
#include "AlleleInfo.h"
#include "utils.h"

using namespace std;

//! A  class to hold VCF info for a single individual
/*!
The constructor parses the line from the VCF and populates the fields. 
*/
class SimpleVCF : public AlleleInfo{
private:
  
    vector<int> countA;
    vector<int> countC;
    vector<int> countG;
    vector<int> countT;
    bool observedGL;
    bool observedPL;
	    
    bool haploidCall;
    //Taken from http://www.broadinstitute.org/gatk/guide/topic?name=intro
    bool unresolvedGT;   //if GT == "./."
    bool homozygousREF;  //if GT == "0/0"
    bool heterozygous;   //if GT == "0/1"
    bool homozygousALT;  //if GT == "1/1"
    bool heterozygousALT;  //if GT == "1/2"
    bool heterozygous2ndALT;  //if GT == "0/2" meaning has ref and second alt
    bool homozygous2ndALT;  //if GT == "0/2" meaning has ref and second alt

    bool deleteCore;

    /* bool resolvedSingleBasePairREF; */
    /* bool resolvedSingleBasePairALT; */
    /* bool allAltResolvedSingleBasePair; */

    //bool closeIndel;
    /* bool isIndel; */
    //Format fields
    string rawFormatNames;
    string rawFormatValues;

    int indexGenotype; //GT
    int indexGenotypeQual; //GQ
    int indexDepth;    //DP
    int indexPL;       //PL
    
    //Format fields
    string formatFieldGT;
    float  formatFieldGQ;
    int    formatFieldDP;
    
    
    string  formatFieldGL;
    string  formatFieldPL;
    int     formatFieldPLHomoRef;
    //for bi-allelic sites
    int     formatFieldPLHetero;
    int     formatFieldPLHomoAlt;
    
    //for tri-allelic sites
    int formatFieldPLHetero1;
    int formatFieldPLHomoAlt1;
    int formatFieldPLHetero2;
    int formatFieldPLHetero12;
    int formatFieldPLHomoAlt2;

    //for penta-allelic sites
    int formatFieldPLHetero3;                        
    int formatFieldPLHetero13;
    int formatFieldPLHetero23;                        
    int formatFieldPLHomoAlt3;
    

    /* unsigned int position; */
    /* string    chrName; */
    /* string    id; */

    /* string    ref; */
    /* string    alt; */
    /* vector<string> altAlleles; */

    /* float    qual;     */
    /* string filter; */

    /* vector<string> fields; */
    /* string infoFieldRaw; */
    /* map<string, string> * infoField; */
    /* bool haveInfoField; */

    vector<string> formatFieldNames;
    vector<string> formatFieldValues;


    //! Returns true if the allele char represented by the char is present, see hasAllele for use with int
    inline bool isThisAllelePresent(char bp) const ;
    CoreVCF *  corevcf;
    void init(const vector<string> & fields, CoreVCF *  corevcf_);
    
public:
    
    //! Constructor 
    /*!
      \param line : The raw line from the VCF file
    */
    SimpleVCF(string line);

    //! Constructor 
    /*!
      \param fields:  The tab delimited fields 
      \param corevcf:  
    */
    SimpleVCF(const vector<string> & fields, CoreVCF *  corevcf_,bool deleteCore_=true); //string line){
    

    //! Dummy constructor, do not use 
    SimpleVCF();

    ~SimpleVCF();

    //! Retrieves a reference to the coreVCF object
    CoreVCF *  getCorevcf();


    //! Retrieves the reference allele as a string as it is in the raw VCF
    string getRef() const;
    //! Retrieves the alternative allele(s) as a string as it is in the raw VCF
    string getAlt() const;
    //! Retrieves the alternative alleles as a vector of strings
    vector<string> getAltAlleles() const;
    //! Retrieves the # of alternative alleles found
    /* int getAltCount() const; */

    /* void parseInfoFields(); */

//! Retrieves the potential associated ID
    string getID() const;
//! Retrieves the chr
    string getChr() const;
//! Retrieves the position on the chr
    unsigned int getPosition() const;
//! Retrieves the filter string 
    string getFilter() const;

//! Retrieves the overall quality
    float  getQual() const;

//! To check if one of the INFO fields exist
/*!
 *   This subroutine is not const because there is a possibility
 *   that it will call parseInfoFields()
 *
  \param tag : Tag of the field in the info line
  \return  : returns true if field exists, false otherwise
  \sa  getInfoField()
*/
    bool hasInfoField(string tag)  ;

    char getOtherAllele(char baseToAvoid) const;
    char getRandomAllele() const;
    char getRandomAlleleUsingPL(int minPLdiffind) const;
    bool isHeterozygousUsingPL(int minPLdiffind) const;


//! To retrieve the information for one of the INFO fields 
/*!
 *
 * The correct invocation has to be used. <int> for int <float> for float etc.. 

  \param tag : Tag of the field in the info line
  \return  : The information associated with the tag, if the tag does not exist, return empty string
  \sa  hasInfoField()
*/
 

    //have we seen the GL field?
    bool getObservedGL() const;
    //have we seen the PL field?
    bool getObservedPL() const;

    //! Retrieves the raw info field, not parsed
    string getInfoFieldRaw() const; 

    //! Retrieves the genotype
    string getGenotype() const;
    //! Retrieves the genotype quality
    float  getGenotypeQual() const;
    //! Retrieves the depth, 
    int    getDepth() const;
    //! Retrieves the depth using the info field (all reads, not just the ones that pass quality)
    int    getDepthInfo() ;


    //! Retrieves the pl field
    string getPL() const;
    //! Retrieves the homo ref likelihood from the pl field
    int    getPLHomoRef() const;
    //! Retrieves the heterozygous likelihood from the pl field
    int    getPLHetero() const;
    //! Retrieves the homo alt likelihood from the pl field
    int    getPLHomoAlt() const;


    //! Sets the close to indel flag
    void    setCloseIndel(bool closeIndel);
    bool    getCloseIndel() const;
    bool    containsIndel() const;

    //! true if the REF allele is A,C,G or T
    bool    isResolvedSingleBasePairREF() const;
    //! true if the ALT allele is A,C,G,T or .
    bool    isResolvedSingleBasePairALT() const;
   //! true if every ALT allele is A,C,G,T or .
    bool    areAllAltResolvedSingleBasePair() const;

    //! true if GT field is "./."
    bool isUnresolvedGT() const; 
    //! true if GT field is "0/0"
    bool isHomozygousREF() const;  
    //! true if GT field is "0/1"
    bool isHeterozygous() const;   
    //! true if GT field is "1/1"
    bool isHomozygousALT() const;
    //! true if GT field is "1/2"
    bool isHeterozygousALT() const;
    //! true if GT field is "0/2"
    bool isHeterozygous2ndALT() const;

    /* friend ostream& operator<<(ostream& os, const SimpleVCF& smvcf); */
    void print(ostream& os) const;

    //! true if the current record contains at least one A
    bool hasAtLeastOneA() const  ;
    //! true if the current record contains at least one C
    bool hasAtLeastOneC() const  ;
    //! true if the current record contains at least one G
    bool hasAtLeastOneG() const  ;
    //! true if the current record contains at least one T
    bool hasAtLeastOneT() const  ;
    
    
    //returns true if this record has this allele (A:1,C:2,G:3,T:4)
    bool hasAllele(int indexAlle) const;

    //! returns the count of reference and alternative allele based on PL field
    pair<int,int> returnLikelyAlleleCountForRefAlt(int minPLdiffind=33) const;
    pair<int,int> returnLikelyAlleleCountForRefAltJustGT() const;

    //returns a string just using GT
    string getAlleCountBasedOnGT() const;

    unsigned int getAltAlleleCount() const;
    int getADforA() const;
    int getADforC() const;
    int getADforG() const;
    int getADforT() const;

    //returns allele count this allele (A:1,C:2,G:3,T:4)
    int getADforAllele(int indexAlle) const;

    template <typename T>
	T   getInfoField(string tag)  {
	/* if(!haveInfoField){ parseInfoFields(); } */
	/* map<string,string>::const_iterator it=infoField->find(tag); */
	/* if(it == infoField->end()){ */
	/*     return T(); */
	/* }else{ */
	return corevcf->getInfoField<T>(tag);
	//}
    }
};
#endif
