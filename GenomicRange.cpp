/*
 * GenomicRange
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GenomicRange.h"

GenomicRange::GenomicRange(string chrName,unsigned int startCoord,unsigned int endCoord){
    if(startCoord>endCoord){
	cerr<<"GenomicRange: Error, the start coordinate  "<<startCoord<<" is greater than the end"<<endl;
	exit(1);
    }

    this->chrName    = chrName;
    this->startCoord = startCoord;
    this->endCoord   = endCoord;

}


GenomicRange::GenomicRange(string chrName,unsigned int startCoord){
    if(startCoord>endCoord){
	cerr<<"GenomicRange: Error, the start coordinate  "<<startCoord<<" is greater than the end"<<endl;
	exit(1);
    }
    this->chrName    = chrName;
    this->startCoord = startCoord;
    this->endCoord   = startCoord;
}



GenomicRange::GenomicRange(){
    this->chrName    = "";
    this->startCoord = 0;
    this->endCoord   = 0;
}

GenomicRange::~GenomicRange(){

}


string GenomicRange::getChrName() const {
    return chrName;
}

void GenomicRange::setChrName(string chrName){
    this->chrName=chrName;
}



unsigned int GenomicRange::getStartCoord() const {
    return startCoord;
}

void GenomicRange::setStartCoord(unsigned int startCoord){
    this->startCoord=startCoord;
}




unsigned int GenomicRange::getEndCoord() const {
    return endCoord;
}

void GenomicRange::setEndCoord(unsigned int endCoord){
    this->endCoord=endCoord;
}




unsigned int GenomicRange::getLength() const {
    return endCoord-startCoord+1;
}

string GenomicRange::asBed() const {
    return ""+chrName+"\t"+stringify(startCoord-1)+"\t"+stringify(endCoord);
}


std::ostream & operator << (std::ostream & os, const GenomicRange & ct){
    os<<ct.chrName<<":"
      <<ct.startCoord<<"-"
      <<ct.endCoord;

    return os;
}
