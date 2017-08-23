/*
 * SingleAllele
 * Date: Mar-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SingleAllele_h
#define SingleAllele_h

#include <iostream>
#include "utils.h"

using namespace std;

class SingleAllele{
    friend class GlacParser;
private:
    int refCount;
    int altCount;
    int totalCount;
    bool isCpg;
public:
    SingleAllele();
    SingleAllele(int refCount,int altCount,bool isCpg);

    SingleAllele(const SingleAllele & other);
    ~SingleAllele();
    /* SingleAllele & operator= (const SingleAllele & other); */

    int getRefCount() const;
    int getAltCount() const;
    bool hasAlt() const;

    int getTotalCount() const;
    bool alleleCountNull() const;

    bool isHeterozygous();
    bool   getIsCpg() const;
    char generateRandomAlleleBias(const char ref,const char alt);
    char generateRandomAllele(const char ref,const char alt);

    int    printEIGENSTRAT();
    char  printPLINK();

    void setRefCount(int refCount);
    void setAltCount(int altCount);
    void addRefCount(int refCount);
    void addAltCount(int altCount);



    void setIsCpg(bool isCpg);
    string toString();

    //SingleAllele &  operator= (const SingleAllele & other);

    SingleAllele & operator= (const SingleAllele & other){
    	refCount   =  other.refCount;
    	altCount   =  other.altCount;
    	totalCount =  other.totalCount;
    	isCpg      =  other.isCpg;
    	return *this;
    }

    friend bool operator== (const SingleAllele & first,const SingleAllele & second);
    friend bool operator!= (const SingleAllele & first,const SingleAllele & second);

    friend ostream& operator<<(ostream& os, const SingleAllele & at){
	os<<at.refCount<<","<<at.altCount<<":"<<at.isCpg;
	return os;
    }


    friend SingleAllele    operator+(const SingleAllele & first,const SingleAllele & second);
    SingleAllele &  operator+=(const SingleAllele & other);

    //friend SingleAllele operator+(const SingleAllele & other);/* { */
    /* 	return SingleAllele( refCount + other.refCount,  */
    /* 			     altCount + other.altCount, */
    /* 			     (isCpg ||  other.isCpg) */
    /* 			     ); */
    /* } */
    
};
#endif
