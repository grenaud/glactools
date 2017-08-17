/*
 * AlleleRecords
 * Date: Mar-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AlleleRecords_h
#define AlleleRecords_h

#include <string>
#include <vector>

#include "SingleAllele.h"
#include "SingleGL.h"

using namespace std;

class AlleleRecords{
private:
    bool glFormat;  

public:
    AlleleRecords(bool glFormat_=false);
    AlleleRecords(uint32_t sizePops_,bool glFormat_);

    AlleleRecords(const AlleleRecords & other);
    ~AlleleRecords();
    AlleleRecords & operator= (const AlleleRecords & other){
	chr        = other.chr;
	chri       = other.chri;
	coordinate = other.coordinate;
	ref        = other.ref;
	alt        = other.alt;
	glFormat   = other.glFormat;

	if(other.glFormat){
	    vectorAlleles = 0;
	    vectorGLs     = new vector<SingleGL>     ( *(other.vectorGLs) );
	}else{
	    vectorAlleles = new vector<SingleAllele> ( *(other.vectorAlleles) );
	    vectorGLs     = 0;
	}
	return *this;
    }
    void copyCoreFields(const AlleleRecords & other);


    string chr;
    uint16_t chri;
    uint32_t coordinate;
    uint32_t sizePops;
    char ref;
    char alt;
    vector<SingleAllele> * vectorAlleles;
    vector<SingleGL>     * vectorGLs;

    //bool buildFromData(char * buffer,
    bool everyRecordNonNull() const;
    bool everyNonChimpAncRecordNonNull() const;
    bool isGlFormat() const;
    void writeBinary(char * buffer) const;
    uint32_t getSizePops() const;

    friend ostream& operator<<(ostream& os, const AlleleRecords & ar){	  
	os<<ar.chr<<"\t";
	os<<stringify(ar.coordinate)<<"\t";
	os<<ar.ref<<",";
	os<<ar.alt<<"\t";
	//exit(1);
	if(ar.glFormat)
	    os<<vectorToString(*(ar.vectorGLs),"\t");
	else
	    os<<vectorToString(*(ar.vectorAlleles),"\t");
	//	exit(1);
	return os;
    }

    friend bool operator== (const AlleleRecords & first,const AlleleRecords & second);

};
#endif
