/*
 * SingleAllele
 * Date: Mar-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SingleAllele.h"

SingleAllele::SingleAllele(){
    refCount=0;
    altCount=0;
    totalCount=0;
    isCpg=false;
}


SingleAllele::SingleAllele(int refCount,int altCount,bool isCpg):
    refCount(refCount),
    altCount(altCount),
    isCpg(isCpg){
    totalCount=refCount+altCount;
}


SingleAllele::SingleAllele(const SingleAllele & other):
    refCount(other.refCount),
    altCount(other.altCount),
    isCpg(other.isCpg){
    totalCount=refCount+altCount;
}

SingleAllele::~SingleAllele(){

}


int SingleAllele::getTotalCount() const{
    return totalCount;
}

bool SingleAllele::alleleCountNull() const{
    return (totalCount==0);
}


int SingleAllele::getRefCount() const{
    return refCount;
}

int SingleAllele::getAltCount() const{
    return altCount;
}

bool SingleAllele::hasAlt() const{
    return (altCount != 0);
}


bool SingleAllele::getIsCpg() const{
    return isCpg;
}


void SingleAllele::setRefCount(int refCount_){
    refCount = refCount_;
    totalCount=refCount+altCount;
}

void SingleAllele::setAltCount(int altCount_){
    altCount = altCount_;
    totalCount=refCount+altCount;
}


void SingleAllele::addRefCount(int refCount_){
    refCount   += refCount_;
    totalCount  = refCount+altCount;
}

void SingleAllele::addAltCount(int altCount_){
    altCount   += altCount_;
    totalCount  = refCount+altCount;
}



void SingleAllele::setIsCpg(bool isCpg_){
    isCpg    = isCpg_;
}

string  SingleAllele::toString(){
    string toReturn =""+ 
	stringify( refCount)+","+
	stringify( altCount)+":"+
	stringify( isCpg);
    return toReturn;
}

SingleAllele  operator+(const SingleAllele & first,const SingleAllele & second){

    return SingleAllele( first.refCount + second.refCount, 
    			 first.altCount + second.altCount,
    			 (first.isCpg   ||  second.isCpg)
    			 );
}


SingleAllele & SingleAllele::operator+=(const SingleAllele & other){
    this->refCount   += other.refCount;
    this->altCount   += other.altCount;
    this->totalCount += other.totalCount;

    this->isCpg    =  (this->isCpg   ||  other.isCpg);
    return *this;
}

int  SingleAllele::printEIGENSTRAT(){
    if(refCount == 2 && altCount == 0) //homo ref
	return 2;

    if(refCount == 1 && altCount == 1) //hetero
	return 1;

    if(refCount == 0 && altCount == 2) //homo alt
	return 0;
    
    if( (refCount == 0 && altCount == 1) ){ //has at least one alt, wrong but we will have to use it
	return 0;	
    }

    if( (refCount == 1 && altCount == 0) ){ //has at least one ref, wrong but we will have to use it
	return 2;
    }
    
    if( (refCount == 0 && altCount == 0) ) //unknown
	return 9;


    //If we reached this position, this means we deal with population allele frequency
    //try to generate something random for populations
    int randIndex1=rand()%totalCount;//returns a number between 0 and (totalCount-1)
    //int randIndex2=rand()%totalCount;//returns a number between 0 and (totalCount-1)
    bool ref1 = (randIndex1<refCount);
    bool ref2 = (randIndex1<refCount);

    if(ref1 && ref2 ) //homo ref
	return 2;
    
    
    if(!ref1 && !ref2 ) //homo alt
	return 0;
    
    return 1;
    
    // cerr<<"SingleAllele.cpp: printEIGENSTRAT() cannot get the EIGENSTRAT for this data record: "<<*this<<" cannot use this for combined individuals (populations)"<<endl;
    // exit(1);
}


bool SingleAllele::isHeterozygous(){
    if(refCount == 0)
	return false;
    if(altCount == 0) 
	return false;
    return true;
}

char  SingleAllele::printPLINK(){
    if(refCount == 2 && altCount == 0) //homo ref
	return 0; //00

    if(refCount == 1 && altCount == 1) //hetero
	return 2; //10

    if(refCount == 0 && altCount == 2) //homo alt
	return 3; //11
    
    if( (refCount == 0 && altCount == 1) ){ //has at least one alt, wrong but we will have to use it
	return 3; //11
    }

    if( (refCount == 1 && altCount == 0) ){ //has at least one ref, wrong but we will have to use it
	return 0; //00
    }
    
    if( (refCount == 0 && altCount == 0) ) //unknown
	return 1; //01


    //If we reached this position, this means we deal with population allele frequency
    //try to generate something random for populations
    int randIndex1=rand()%totalCount;//returns a number between 0 and (totalCount-1)
    //int randIndex2=rand()%totalCount;//returns a number between 0 and (totalCount-1)
    bool ref1 = (randIndex1<refCount);
    bool ref2 = (randIndex1<refCount);

    if(ref1 && ref2 ) //homo ref
	return 0; //00
    
    
    if(!ref1 && !ref2 ) //homo alt
	return 3; //11
    
    return 1;     //01
}

void SingleAllele::downsample(const int c){
    if(totalCount==0)
	return ;
    int refCount_   = 0;
    int altCount_   = 0;
    int totalCount_ = 0;
    
    
    for(int i=0;i<c;i++){
	int r=randomInt(1,totalCount);
	if(r<=refCount){//sampled ref
	    refCount_++;
	    refCount--;	    
	}else{ //sampled alt
	    altCount_++;
	    altCount--;	    	    
	}
    }
    refCount=refCount_;
    altCount=altCount_;
    totalCount=(refCount+altCount);    
}

char SingleAllele::generateRandomAlleleBias(const char ref,const char alt){
    if(totalCount==0){
	// cerr<<"SingleAllele: cannot call generateRandomAlleleBias() where the allele count is 0 :"<<*this<<endl;
	// exit(1);
	return 'N';
    }
    int randIndex=rand()%totalCount;//returns a number between 0 and (totalCount-1)

    if(randIndex<refCount){
	return ref;
    }else{
	return alt;
    }

}

char SingleAllele::generateRandomAllele(const char ref,const char alt){
    if(totalCount==0){
	// cerr<<"SingleAllele: cannot call generateRandomAlleleBias() where the allele count is 0 :"<<*this<<endl;
	// exit(1);
	return 'N';
    }

    if(refCount > 0  && altCount == 0) //ref
	return ref;

    if(refCount == 0 && altCount > 0) //alt
	return alt;

    if(randomBool())
	return ref;
    else
    	return alt;
}


bool operator== (const SingleAllele & first,const SingleAllele & second){
    return (  (first.refCount   ==  second.refCount) &&
	      (first.altCount   ==  second.altCount) &&
	      (first.totalCount ==  second.totalCount) &&
	      (first.isCpg      ==  second.isCpg)  );	      
}


bool operator!= (const SingleAllele & first,const SingleAllele & second){
    return !(first==second);
}


// SingleAllele & SingleAllele::operator= (const SingleAllele & other){
//     refCount   =  other.refCount;
//     altCount   =  other.altCount;
//     totalCount =  other.totalCount;
//     isCpg      =  other.isCpg;
//     return *this;
// }
