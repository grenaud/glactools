/*
 * CoreVCF
 * Date: Sep-30-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "CoreVCF.h"

CoreVCF::CoreVCF(const vector<string> & fields){
    closeIndel=false;

    fieldIndex=0;

    // if(fields.size() == 11){

    // 	if( !fields[10].empty()  ){
    // 	    cerr<<"CoreVCF ERROR, found # fields =  ("<<fields.size()<<") for : #"<<vectorToString(fields,"\t")<<"#"<<endl;
    // 	    for(unsigned int i=0;i<fields.size();i++){
    // 		cerr<<"fields["<<i<<"] = "<<fields[i]<<endl;
    // 	    }
    // 	    exit(1);
    // 	}
	
    // }else{

    if(fields.size() < 9){
	cerr<<"CoreVCF ERROR, not enough fields ("<<fields.size()<<") for : #"<<vectorToString(fields,"\t")<<"#"<<endl;
	for(unsigned int i=0;i<fields.size();i++){
	    cerr<<"fields["<<i<<"] = "<<fields[i]<<endl;
	}
	exit(1);
    }

	//    }

    chrName=                     fields[fieldIndex++];
    position=string2uint(        fields[fieldIndex++]);

    if(fields.size() != 9 )
    	id     =                 fields[fieldIndex++];
    else
    	id     =                 "NA";


    ref    =                     fields[fieldIndex++];
    alt    =                     fields[fieldIndex++];

    //for sites with multiple alt bases
    altAlleles      =     allTokens(alt,',');
    allAltResolvedSingleBasePair=true;
    for(unsigned int i=0;i<altAlleles.size();i++){
	allAltResolvedSingleBasePair&=validAltBP( altAlleles[i] ); //true if ref = A,C,G,T or .
    }
    
    
    //boolean flags for insert
    isIndel=isInsert(ref) || isInsert(alt);
    
    //cerr<<"isIndel "<<isIndel<<endl;
    //boolean flags for a single bp in ref or alt
    resolvedSingleBasePairREF=validOneBP(ref); //true if ref = A,C,G or T
    resolvedSingleBasePairALT=validAltBP(alt); //true if ref = A,C,G,T or .
    
    int qualFieldIdx=fieldIndex++;
    if(fields[qualFieldIdx] == "."){
	qual=0.0;
    }else{
	qual   = destringify<float>( fields[qualFieldIdx] );
    }
    filter =                     fields[fieldIndex++];


    //INFO FIELD
    fieldIndexINFO    = fieldIndex;
    infoFieldRaw      = fields[fieldIndex++] ;
    haveInfoField=false;
    //if INDEL MARKED
    isIndel = isIndel || strBeginsWith(infoFieldRaw,"INDEL");


    //increasing to skip GT FIELD
    fieldIndex++;

}

CoreVCF::~CoreVCF(){
    // cerr<<"CoreVCF des1"<<endl;
    if(haveInfoField)
	delete infoField;
    // cerr<<"CoreVCF des2"<<endl;
}



void CoreVCF::parseInfoFields(){
    haveInfoField=true;
    infoField    =        info2map( infoFieldRaw );
}

bool   CoreVCF::hasInfoField(string tag) {
    if(!haveInfoField){ parseInfoFields(); }
    return (infoField->find(tag)  != infoField->end());
}


int     CoreVCF::getDepthInfo() {
    if(!haveInfoField){ parseInfoFields(); }

    if(hasInfoField("DP")){
	int toReturn = getInfoField<int>("DP");
	return toReturn;
    }
    return -1;
}



string CoreVCF::getRef() const{
    return ref;
}

string CoreVCF::getAlt() const{
    return alt;
}

string CoreVCF::getID() const{
    return id;
}

string CoreVCF::getChr() const{
    return chrName;
}

unsigned int CoreVCF::getPosition() const{
    return position;
}

string  CoreVCF::getFilter() const{
    return filter;
}

string CoreVCF::getInfoFieldRaw() const{
    return infoFieldRaw;
}

int CoreVCF::getFieldIndexAndIncrease(){    
    return fieldIndex++;
}


int CoreVCF::getFieldIndexINFO(){
    return fieldIndexINFO;
}



float     CoreVCF::getQual() const{
    return qual;
}

void    CoreVCF::setCloseIndel(bool closeIndel){
    this->closeIndel=closeIndel;
}

bool CoreVCF::getCloseIndel() const{
    return closeIndel;
}

bool CoreVCF::containsIndel() const{
    return isIndel;
}

bool    CoreVCF::isResolvedSingleBasePairREF() const{
    return resolvedSingleBasePairREF;
}

bool    CoreVCF::isResolvedSingleBasePairALT() const{
    return resolvedSingleBasePairALT;
}
