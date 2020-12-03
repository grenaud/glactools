/*
 * test
 * Date: Jul-31-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include "VCFreader.h"
#include "SimpleVCF.h"

#include "libgab.h"

using namespace std;

int main (int argc, char *argv[]) {

    
    string filename = string(argv[1]);
    VCFreader     vcfr (filename,5);



    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();
	cout<<toprint->getChr()<<"\t"<<toprint->getPosition()<<"\t"<<toprint->getRef()<<"\t"<<toprint->getAlt()<<"\t"<<toprint->isCpg()<<endl;
    }

    return 0;
}

