/*
 * readACFFile
 * Date: Sep-27-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 * 
 * 
 */

#include <iostream>
#include <fstream>

#include "GlacParser.h"
#include "AlleleRecords.h"


using namespace std;

int main (int argc, char *argv[]) {


    string glacfile  = string(argv[argc-1]);
    cerr<<"reading: "<<glacfile<<endl;

    GlacParser gp (glacfile);

    AlleleRecords * arr;

    while(gp.hasData()){	
	arr = gp.getData();	
	cout<<*arr<<endl;
    }

    cerr<<"done!"<<endl;
    
    return 0;
}

