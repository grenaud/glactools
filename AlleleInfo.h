/*
 * AlleleInfo
 * Date: Aug-29-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef AlleleInfo_h
#define AlleleInfo_h

#include <iostream>
#include <string>

using namespace std;

//! A  virtual class for a classes that hold allele data 
/*!
This class contains basic subroutines for the retrieval of 
allele data.
*/
class AlleleInfo{

 private:
    bool isCpg_;
 protected: //derived classes can have access
    int typeOfData; //1 = vcf, 2 = bamtable
 public:
    AlleleInfo(){
	isCpg_=false;    
    }
    virtual ~AlleleInfo() {   //cout<<"DESTRUCTOR AlleleInfo"<<endl;     
    };
    virtual unsigned int getPosition() const =0;
    virtual string       getChr() const =0;
    virtual char         getRandomAllele() const =0;
    /* virtual ostream& operator<<(ostream& os) =0; */

    virtual void print(ostream& os) const = 0 ;

    virtual bool hasAtLeastOneA() const = 0 ;
    virtual bool hasAtLeastOneC() const = 0 ;
    virtual bool hasAtLeastOneG() const = 0 ;
    virtual bool hasAtLeastOneT() const = 0 ;
    virtual bool hasAllele(int indexAlle)      const = 0 ;


    /* ostream& operator<< (ostream& out, const AlleleInfo& al) ; */
    void setCpg(bool flag){
	isCpg_=flag;
    }

    bool isCpg() const{
	return isCpg_;
    }
    


    friend ostream& operator<<(ostream& str, AlleleInfo const& al){
	al.print(str);
	return str;
    }

};
#endif
