/*
 * SingleGL
 * Date: Mar-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SingleGL_h
#define SingleGL_h

#include <iostream>
#include "utils.h"

using namespace std;

class SingleGL{
private:

    uint8_t rrGL;
    uint8_t raGL;
    uint8_t aaGL;
    bool isCpg;
public:
    SingleGL();
    SingleGL(uint8_t rrGL,uint8_t raGL,uint8_t aaGL,bool isCpg);

    SingleGL(const SingleGL & other);
    ~SingleGL();
    /* SingleGL & operator= (const SingleGL & other); */

    uint8_t getrrGL();
    uint8_t getraGL();
    uint8_t getaaGL();

    /* int getTotalCount(); */
    /* bool alleleCountNull(); */

    /* bool isHeterozygous(); */
    bool   getIsCpg();
    /* char generateRandomAlleleBias(const char ref,const char alt); */
    /* char generateRandomAllele(const char ref,const char alt); */

    /* int    printEIGENSTRAT(); */
    /* char  printPLINK(); */

    /* void setRefCount(int refCount); */
    /* void setAltCount(int altCount); */
    /* void addRefCount(int refCount); */
    /* void addAltCount(int altCount); */



    void setIsCpg(bool isCpg);
    string toString();

    //SingleGL &  operator= (const SingleGL & other);

    /* SingleGL & operator= (const SingleGL & other){ */
    /* 	refCount   =  other.refCount; */
    /* 	altCount   =  other.altCount; */
    /* 	totalCount =  other.totalCount; */
    /* 	isCpg      =  other.isCpg; */
    /* 	return *this; */
    /* } */

    friend bool operator== (const SingleGL & first,const SingleGL & second);
    friend bool operator!= (const SingleGL & first,const SingleGL & second);

    friend ostream& operator<<(ostream& os, const SingleGL & at){
	//os<<"op<"<<endl;
	/* cerr<<at.rrGL<<","<<at.raGL<<","<<at.aaGL<<":"<<at.isCpg; */
	/* exit(1); */
	//os<<at.rrGL<<","<<at.raGL<<","<<at.aaGL<<":"<<at.isCpg;
	return os;
    }


    friend SingleGL    operator+(const SingleGL & first,const SingleGL & second);
    SingleGL &  operator+=(const SingleGL & other);

    
};
#endif
