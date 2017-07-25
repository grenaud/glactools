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

    AlleleRecords(const AlleleRecords & other);
    ~AlleleRecords();
    AlleleRecords & operator= (const AlleleRecords & other){
	chr        = other.chr;
	coordinate = other.coordinate;
	ref        = other.ref;
	alt        = other.alt;
	vectorAlleles = new vector<SingleAllele> ( *(other.vectorAlleles) );
	return *this;
    }
    

    string chr;
    uint16_t chri;
    unsigned int coordinate;
    char ref;
    char alt;
    vector<SingleAllele> * vectorAlleles;
    vector<SingleGL>     * vectorGLs;

    //bool buildFromData(char * buffer,
    bool everyRecordNonNull() const;
    bool everyNonChimpAncRecordNonNull() const;


    friend ostream& operator<<(ostream& os, const AlleleRecords & ar){	  
	os<<ar.chr<<"\t";
	os<<stringify(ar.coordinate)<<"\t";
	os<<ar.ref<<",";
	os<<ar.alt<<"\t";

	if(ar.glFormat)
	    os<<"a"<<vectorToString(*(ar.vectorGLs),"\t");
	else
	    os<<vectorToString(*(ar.vectorAlleles),"\t");
	//	exit(1);
	return os;
    }

    friend bool operator== (const AlleleRecords & first,const AlleleRecords & second);

};
#endif
