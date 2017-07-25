/*
 * GenomicRange
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#ifndef GenomicRange_h
#define GenomicRange_h

#include <iostream>
#include <string>
#include <stdlib.h>

#include "utils.h"

using namespace std;

class GenomicRange{

private:
    string chrName;
    unsigned int startCoord;
    unsigned int endCoord;

public:
    GenomicRange();
    GenomicRange(string chrName,unsigned int startCoord);
    GenomicRange(string chrName,unsigned int startCoord,unsigned int endCoord);
    ~GenomicRange();
	
    string getChrName() const;
    unsigned int getStartCoord() const;
    unsigned int getEndCoord() const;

    void setChrName(string chrName);
    void setStartCoord(unsigned int startCoord);
    void setEndCoord(unsigned int endCoord);
    string asBed() const ;

    unsigned int getLength() const;


    friend std::ostream & operator<<(std::ostream & os, const GenomicRange & ct);



};
#endif
