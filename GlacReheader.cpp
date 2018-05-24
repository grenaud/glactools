/*
 * GlacReheader
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacReheader.h"

GlacReheader::GlacReheader(){

}

GlacReheader::~GlacReheader(){

}


string GlacReheader::usage() const{

    
    return string("") +"glactools reheader [options] <glac file> <text file>\n"+
	"This program will remove the header from the <glac file> with <text file> and will print to STDOUT\n"+
	"WARNING: Make sure the header is valid as no error checking is done\n"+
	"         and downstream programs might not work\n"+	
	"ex:  glactools reheader data.acf.gz header.txt"+"\n"+
	"Options:"+"\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\n"
	;
}



int GlacReheader::run(int argc, char *argv[]){

    int lastOpt=1;

    for(int i=1;i<(argc-1);i++){ 
        //cout<<i<<"\t"<<string(argv[i])<<endl;
        if((string(argv[i]) == "-")  ){
            lastOpt=i;
            break;          
        }
        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }

        if(string(argv[i]) == "-u"){
            uncompressed=true;
            continue;
        }


        cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    if(lastOpt != (argc-2)){
        cerr<<"The last arguments are the  <glf|acf file> <text file> "<<endl;
        return 1;               
    }

    
    string fileglac                    = string(argv[lastOpt  ]);
    string newheader                   = string(argv[lastOpt+1]);
    GlacParser gp (fileglac);

    igzstream myFile;
    myFile.open(newheader.c_str(), ios::in);
    stringstream header;
    string line;
    if (myFile.good()){
        while ( getline (myFile,line)){   
	    header<<line<<endl;
	}
	myFile.close();
    }else{
        cerr << "GlacReheader: Unable to open file "<<newheader<<endl;
        return 1;
    }
    

    GlacWriter * gw = new GlacWriter(gp.getSizePops(),
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     1,//compression threads
				     uncompressed);

    
    if(!gw->writeHeader(header.str())){
	cerr<<"GlacRemovePop: error writing header "<<endl;
	exit(1);
    }
	

    AlleleRecords * arr;

    while(gp.hasData()){
	arr = gp.getData();

	if(!gw->writeAlleleRecord(arr)){
     	    cerr<<"GlacReheader: error record "<<*arr<<endl;
     	    exit(1);
     	}
    }


    delete(gw);



    return 0;
}

