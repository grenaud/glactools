/*
 * CoreVCF
 * Date: Sep-30-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef CoreVCF_h
#define CoreVCF_h

#include <iostream>
#include <string>
#include <vector> 
#include <map>
#include "utils.h"

using namespace std;

class CoreVCF{
 private:
    bool closeIndel;
    unsigned int position;
    string    chrName;
    string    id;

    string    ref;
    string    alt;
    vector<string> altAlleles;
    bool allAltResolvedSingleBasePair;
    bool isIndel;

    bool resolvedSingleBasePairREF;
    bool resolvedSingleBasePairALT;
    float    qual;    
    string filter;


    vector<string> fields;
    string infoFieldRaw;    
    map<string, string> * infoField;
    bool haveInfoField;
    int fieldIndex;
    int fieldIndexINFO;
    string rawFormatNames;
    vector<string> * formatNames;

    template <typename T>
	T   getInfoField(string tag)  {
	if(!haveInfoField){ parseInfoFields(); }
	map<string,string>::const_iterator it=infoField->find(tag);
	if(it == infoField->end()){
	    return T();
	}else{
	    return destringify<T>(it->second);
	}
    }


    void parseInfoFields();


    bool   hasInfoField(string tag) ;    
    int     getDepthInfo() ;   
    string getRef() const;
    string getAlt() const;
    string getID() const;
    string getChr() const;
    unsigned int getPosition() const;
    string  getFilter() const;
    string getInfoFieldRaw() const;

    float     getQual() const;
    void    setCloseIndel(bool closeIndel);
    bool getCloseIndel() const;
    bool containsIndel() const;
    bool    isResolvedSingleBasePairREF() const;
    bool    isResolvedSingleBasePairALT() const;

    int getFieldIndexAndIncrease();
    int getFieldIndexINFO();

public:

    CoreVCF(const vector<string> & fields);
    CoreVCF(const CoreVCF & other);
    ~CoreVCF();
    CoreVCF & operator= (const CoreVCF & other);
    const vector<string> * getFormatNames();
    friend class SimpleVCF;
};
#endif
