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




uint8_t SingleGL::getrrGL(){
    return rrGL;
}

uint8_t SingleGL::getraGL(){
    return raGL;
}

uint8_t SingleGL::getaaGL(){
    return aaGL;
}

bool SingleGL::getIsCpg(){
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
