/*
 * AlleleRecords
 * Date: Mar-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "AlleleRecords.h"

AlleleRecords::AlleleRecords(bool glFormat_){
    chr="NA";
    coordinate=0;
    ref='N';
    alt='N';
    vectorAlleles=0;
    vectorGLs=0;
    glFormat = glFormat_;
}

AlleleRecords::AlleleRecords(const AlleleRecords & other){
    chr        = other.chr;
    coordinate = other.coordinate;
    ref        = other.ref;
    alt        = other.alt;
    vectorAlleles = new vector<SingleAllele> ( *(other.vectorAlleles) );
    vectorGLs     = new vector<SingleGL>     ( *(other.vectorGLs) );

}

AlleleRecords::~AlleleRecords(){
    if(vectorAlleles != 0)
	delete(vectorAlleles);
    if(vectorGLs != 0)
	delete(vectorGLs);
}


bool AlleleRecords::everyRecordNonNull() const{
    bool toReturn=false;
    for(unsigned int i=0;i<vectorAlleles->size();i++){
	toReturn= toReturn || (vectorAlleles->at(i).getTotalCount() == 0 );
    }
    return !toReturn;
}


bool AlleleRecords::everyNonChimpAncRecordNonNull() const{
    bool toReturn=false;
    for(unsigned int i=2;i<vectorAlleles->size();i++){
	toReturn= toReturn || (vectorAlleles->at(i).getTotalCount() == 0 );
    }
    return !toReturn;
}



bool operator== (const AlleleRecords & first,const AlleleRecords & second){
    bool normalData= (  (first.chr        ==  second.chr) &&
			(first.ref        ==  second.ref) &&
			(first.alt        ==  second.alt)   );	      
    if(!normalData)
	return false;

    if(first.vectorAlleles->size() != second.vectorAlleles->size()){
	return false;
    }

    if(first.glFormat != second.glFormat)
	return false;

    if(first.glFormat){
	for(unsigned int i=0;i<first.vectorAlleles->size();i++){
	    if(first.vectorAlleles->at(i) !=  second.vectorAlleles->at(i) ){
		return false;
	    }
	}
    }else{
	for(unsigned int i=0;i<first.vectorGLs->size();i++){
	    if(first.vectorGLs->at(i) !=  second.vectorGLs->at(i) ){
		return false;
	    }
	}
    }

    return true;
}
