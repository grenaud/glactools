/*
 * GlacViewer
 * Date: Jul-25-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacViewer.h"
#include "GlacParser.h"

GlacViewer::GlacViewer(){

}

GlacViewer::~GlacViewer(){

}



string GlacViewer::usage() const{
        
    return string(string("") +"view <options> [gl|ac file]");

}


int GlacViewer::run(int argc, char *argv[]){
    string glacfile  = string(argv[argc-1]);
    GlacParser gp (glacfile);
    AlleleRecords * test;

    while(gp.hasData()){

	test = gp.getData();
	cout<<*test<<endl;
	//cout<<test<<endl;
    }

    return 0;
}

