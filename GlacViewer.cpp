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
    int lastOpt=1;

    //for(int i=1;i<(argc);i++){ 

    

    for(int i=1;i<(argc);i++){ 
	//cout<<i<<"\t"<<string(argv[i])<<endl;

        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }

        // if(string(argv[i]) == "-u"){
        //     uncompressed=true;
        //     continue;
        // }
	
	if(string(argv[i]) == "-h"){
	    printheader=true;
            continue;
        }

        cerr<<"Error: unknown option "<<string(argv[i])<<endl;
        return 1;	
    }
    
    string glacfile  = string(argv[argc-1]);
    GlacParser gp (glacfile);
    AlleleRecords * test;
    if(printheader)
	cout<<gp.getHeader()<<endl;
    while(gp.hasData()){

	test = gp.getData();
	// cout<<"run"<<endl;
	// cout<<test<<endl;
	cout<<*test<<endl;

    }

    return 0;
}

