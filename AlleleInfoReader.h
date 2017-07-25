/*
 * AlleleInfoReader
 * Date: Aug-29-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AlleleInfoReader_h
#define AlleleInfoReader_h

#include <iostream>
#include <memory>

#include "AlleleInfo.h"

using namespace std;


//! A  virtual class for a classes that read allele data 
/*!
  This class contains basic subroutines for classes that read allele data (VCF, BAMtable)
*/
class AlleleInfoReader{
private:

public:
    
    virtual ~AlleleInfoReader() {    //cout<<"DESTRUCTOR Allele"<<endl; 
    };
    virtual bool hasData() =0 ;
    virtual AlleleInfo * getData() =0;
    /* virtual auto_ptr<AlleleInfo>  getData() =0; */
    virtual  void repositionIterator(string chrName,int start,int end) =0;
    
    
};
#endif
