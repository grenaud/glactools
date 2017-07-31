/*
 * SingleGL
 * Date: Mar-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SingleGL.h"

SingleGL::SingleGL(){
    rrGL=0;
    raGL=0;
    aaGL=0;
    isCpg=false;
}


SingleGL::SingleGL(uint8_t rrGL,uint8_t raGL,uint8_t aaGL,bool isCpg):
    rrGL(rrGL),
    raGL(raGL),
    aaGL(aaGL),
    isCpg(isCpg){

}


SingleGL::SingleGL(const SingleGL & other):
    rrGL(other.rrGL),
    raGL(other.raGL),
    aaGL(other.aaGL),
    isCpg(other.isCpg){
}

SingleGL::~SingleGL(){

}




uint8_t SingleGL::getrrGL() const{
    return rrGL;
}

uint8_t SingleGL::getraGL() const{
    return raGL;
}

uint8_t SingleGL::getaaGL() const{
    return aaGL;
}

void SingleGL::setrrGL(uint8_t rrGL_){
    rrGL = rrGL_;
}

void SingleGL::setraGL(uint8_t raGL_){
    raGL = raGL_;
}

void SingleGL::setaaGL(uint8_t aaGL_){
    aaGL = aaGL_;
}

bool SingleGL::getIsCpg() const{
    return isCpg;
}




void SingleGL::setIsCpg(bool isCpg_){
    isCpg    = isCpg_;
}

string  SingleGL::toString(){
    string toReturn =""+ 
	stringify( rrGL)+","+
	stringify( raGL)+","+
	stringify( aaGL)+":"+
	stringify( isCpg);
    return toReturn;
}

bool SingleGL::alleleCountNull() const{
    if ( (rrGL == 0) &&
	 (raGL == 0) &&
	 (aaGL == 0) )
	return true;

    if ( (rrGL == UINT8_MAX) &&
	 (raGL == UINT8_MAX) &&
	 (aaGL == UINT8_MAX) )
	return true;
	 
    return false;
}

pair<int,int> SingleGL::returnLikelyAlleleCountForRefAlt(int minPLdiffind) const{

        
    if ( (raGL-rrGL) >= minPLdiffind && (aaGL-rrGL) >= minPLdiffind) {  //high likelihood of homo ref, produce 2 alleles ref
	// refAlleles+=2;
	// altAlleles+=0;
	return pair<int,int>(2,0);
    } else{
	if ((raGL-aaGL) >= minPLdiffind && (rrGL-aaGL) >= minPLdiffind) {  //high likelihood of homo alt, produce 2 alleles alt
	    // refAlleles+=0;
	    // altAlleles+=2;
	    return pair<int,int>(0,2);
	} else {
	    if ((rrGL-raGL) >= minPLdiffind && (aaGL-raGL) >= minPLdiffind) { //high likelihood of hetero, produce 1 allele of each
		// refAlleles+=1;
		// altAlleles+=1;
		return pair<int,int>(1,1);
	    }else{
		if ((rrGL-aaGL) >= minPLdiffind && (raGL-aaGL) <minPLdiffind ) { //high likelihood of at least one alt, produce 1 allele alt 
		    // refAlleles+=0;
		    // altAlleles+=1;
		    return pair<int,int>(0,1);
		}else{
		    if ( (aaGL-rrGL) >= minPLdiffind && (raGL-rrGL) < minPLdiffind ) { // high likelihood of at least one ref, produce 1 allele ref
			// refAlleles+=1;
			// altAlleles+=0;
			return pair<int,int>(1,0);
		    }else{
			// refAlleles+=0;
			// altAlleles+=0;
			return pair<int,int>(0,0);
		    }
		}
	    }
	}
    }// end all cases

    return pair<int,int>(0,0);
}



bool operator== (const SingleGL & first,const SingleGL & second){
    return (  (first.rrGL   ==  second.rrGL) &&
	      (first.raGL   ==  second.raGL) &&	  
	      (first.aaGL   ==  second.aaGL) &&	      
	      (first.isCpg      ==  second.isCpg)  );	      
}


bool operator!= (const SingleGL & first,const SingleGL & second){
    return !(first==second);
}

// SingleGL & SingleGL::operator= (const SingleGL & other){
//     refCount   =  other.refCount;
//     altCount   =  other.altCount;
//     totalCount =  other.totalCount;
//     isCpg      =  other.isCpg;
//     return *this;
// }
