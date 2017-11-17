/*
 * CoreVCF
 * Date: Sep-30-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "CoreVCF.h"

CoreVCF::CoreVCF(){
    closeIndel=false;
    fieldIndex=0;
}

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

    //chrName=                     fields[fieldIndex++];
    setName( fields[fieldIndex++].c_str() );

    //position=string2uint(        fields[fieldIndex++]);
    setPos( fields[fieldIndex++].c_str() );

    //if(fields.size() != 9 )
    //id     =                 fields[fieldIndex++];
    setID( fields[fieldIndex++].c_str() );

    // else
    // 	id     =                 "NA";

    //cerr<<"CoreVCF id="<<id<<" "<<fieldIndex<<endl;
    // ref    =                     fields[fieldIndex++];
    // resolvedSingleBasePairREF=validOneBP(ref); //true if ref = A,C,G or T
    setREF( fields[fieldIndex++].c_str() );


    setALT( fields[fieldIndex++].c_str() );
    setQUAL(  fields[fieldIndex++].c_str() );
    // int qualFieldIdx=fieldIndex++;
    // //cerr<<"CoreVCF qualFieldIdx="<<qualFieldIdx<<" "<<fieldIndex<<endl;
    // if(fields[qualFieldIdx] == "."){
    // 	qual=0.0;
    // }else{
    // 	qual   = destringify<float>( fields[qualFieldIdx] );
    // }
    //cerr<<"CoreVCF qual="<<qual<<endl;

    //filter =                     fields[fieldIndex++];
    setFILTER(  fields[fieldIndex++].c_str() );
    //INFO FIELD
    fieldIndexINFO    = fieldIndex;
    setINFO(  fields[fieldIndex++].c_str() );
    // cerr<<"CoreVCF fieldIndex "<<fieldIndex<<endl;    
    // cerr<<"CoreVCF filter="<<filter<<endl;


    //cerr<<"CoreVCF fieldIndex0 "<<fieldIndex<<endl;    


    //infoFieldRaw      = fields[fieldIndex++] ;
    // cerr<<"CoreVCF info="<<infoFieldRaw<<endl;
    // cerr<<"CoreVCF fieldIndex1 "<<fieldIndex<<endl;    
    //infoFieldsNames   = allTokens(infoFieldRaw ,':');



    //cerr<<"CoreVCF fieldIndex2 "<<fieldIndex<<endl;    
    //increasing to skip GT FIELD
    
    setFORMAT(  fields[fieldIndex++].c_str() );


    //cerr<<"CoreVCF format="<<rawFormatNames<<endl;
    //    fieldIndex++;

}

CoreVCF::~CoreVCF(){
    // cerr<<"CoreVCF des1"<<endl;
    if(haveInfoField)
	delete infoField;
    // cerr<<"CoreVCF des2"<<endl;
    delete formatNames;
}

void CoreVCF::setName(  const char * p){
    chrName=                 string(p);
}

void CoreVCF::setPos(   const char * p){
    position = (unsigned int)strtoul(p, NULL, 0);
    //position=string2uint(        fields[fieldIndex++]);
}

void CoreVCF::setID(    const char * p){
    id     =                 string(p);
}

void CoreVCF::setREF(   const char * p){
    ref                      = string(p);
    resolvedSingleBasePairREF= validOneBP(ref); //true if ref = A,C,G or T
}

void CoreVCF::setALT(   const char * p){
    alt    =                     string(p);
    // cerr<<"CoreVCF ref="<<ref<<" "<<fieldIndex<<endl;
    // cerr<<"CoreVCF alt="<<alt<<" "<<fieldIndex<<endl;
    // cerr<<"CoreVCF fieldIndex "<<fieldIndex<<endl;    
    //for sites with multiple alt bases
    altAlleles      =     allTokens(alt,',');
    //cerr<<"CoreVCF fieldIndex "<<fieldIndex<<endl;    
    allAltResolvedSingleBasePair=true;
    //cerr<<"CoreVCF fieldIndex "<<fieldIndex<<endl;    
    for(unsigned int i=0;i<altAlleles.size();i++){
	allAltResolvedSingleBasePair&=validAltBP( altAlleles[i] ); //true if ref = A,C,G,T or .
    }
    //cerr<<"CoreVCF fieldIndex "<<fieldIndex<<endl;        
    
    //boolean flags for insert
    isIndel=isInsert(ref) || isInsert(alt);
    
    //cerr<<"isIndel "<<isIndel<<endl;
    //boolean flags for a single bp in ref or alt
    resolvedSingleBasePairALT=validAltBP(alt); //true if ref = A,C,G,T or .
    //cerr<<"CoreVCF fieldIndex "<<fieldIndex<<endl;    

}


void CoreVCF::setQUAL(  const char * p){
    //int qualFieldIdx=fieldIndex++;
    //cerr<<"CoreVCF qualFieldIdx="<<qualFieldIdx<<" "<<fieldIndex<<endl;
    //if(fields[qualFieldIdx] == "."){
    if(strcmp(p,".") == 0 ){
	qual=0.0;
    }else{
	qual   = atof(p); //destringify<float>( fields[qualFieldIdx] );
    }

}

void CoreVCF::setFILTER(const char * p){
    filter =           string(p);
}

void CoreVCF::setINFO(  const char * p){

    infoFieldRaw      = string(p);
    haveInfoField=false;
    //if INDEL MARKED
    isIndel = isIndel || strBeginsWith(infoFieldRaw,"INDEL");

}

void CoreVCF::setFORMAT(const char * p){
    rawFormatNames  =                string(p);
    formatNames     = new vector<string> (allTokens(rawFormatNames ,':'));
}


const vector<string> * CoreVCF::getFormatNames(){
    return formatNames;
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
