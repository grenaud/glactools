/*
 * GlacCAT
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacCAT.h"

GlacCAT::GlacCAT(){

}

GlacCAT::~GlacCAT(){

}


string GlacCAT::usage() const{

    
    return string(string("glactools") +" cat [options] <glf file1> <glf file2> .. "+"\n"+
                  "\nThis program concatenates ACF/GLF files given that they were from the same genome assembly\n"+    
		  "uses the first file as header and prints to STDOUT\n"
		  "\n"+
		  "Options:\n"+
		  "\t"+"-u" + "\t\t\t"+"Unsafe mode, will not check if the acf/glf files are compatible\n"+
		  "\t"+"  " + "\t\t\t"+"use this if you are sure the acf come from the same reference genome, it is much faster (default: "+booleanAsString(false)+")\n"
		  );
}



int GlacCAT::run(int argc, char *argv[]){


    AlleleRecords * arr;
    //AlleleRecords * arw;
    GlacWriter * gw=NULL;
    string sqLines="";
    string defline="";
    uint16_t chriRead=0;
    bool unsafemode=false;
    int firstIndex=1;
    if(string(argv[1]) == "-u"){
	unsafemode=true;
	firstIndex=2;
    }

    for(int i=firstIndex;i<(argc);i++){ 
	
	string filename = string(argv[i]);
	GlacParser gp (filename);
	size_t sizeRecords;

	if(i==firstIndex){//first file
	    gw = new GlacWriter(gp.getSizePops(),
				false,
				2,
				1,//compression threads
				uncompressed);

	    string newheader=gp.getHeader();
	    sqLines = gp.getHeaderSQ();
	    defline = gp.getDefline();

	    if(!gw->writeHeader(newheader)){
		cerr<<"GlacCat: error writing header "<<endl;
		return 1;
	    }
	    
	    //first record
	    if(gp.hasData()){
		arr = gp.getData();
		chriRead = arr->chri;
		if(!gw->writeAlleleRecord(arr)){
		    cerr<<"GlacCAT: error writing record "<<arr<<endl;
		    exit(1);
		}
	    }


	    if(unsafemode){
		sizeRecords=gp.getSizeRecord();
		char * buffer;
		while(true){
		    buffer = gp.fillBuffer(sizeRecords);
		    if(buffer == 0){
			break;
		    }
		    if(!gw->writeBuffer(buffer,sizeRecords)){
			cerr<<"GlacCAT: error writing record "<<endl;
			exit(1);
		    }
		    delete(buffer);
		}
	    }else{	    
		while(gp.hasData()){
		    arr = gp.getData();

		    if(!gw->writeAlleleRecord(arr)){
			cerr<<"GlacCAT: error writing record "<<arr<<endl;
			exit(1);
		    }
		}
	    }
	}else{
	    if(sqLines != gp.getHeaderSQ()){
		cerr<<"GlacCat: error the SQ lines do not match in headers in file: "<<string(argv[i])<<" does not match the ones in "<<string(argv[1])<<endl;
		return 1;
	    }
	    if(defline != gp.getDefline()){
		cerr<<"GlacCat: error the defline do not match in headers in file: "<<string(argv[i])<<" does not match the one in "<<string(argv[1])<<endl;
		return 1;
	    }

	    //first record
	    if(gp.hasData()){
		arr = gp.getData();
		if(chriRead > arr->chri){
		    cerr<<"GlacCAT: WARNING chromosomal coordinate in  "<<arr<<" is greater than in the previous file, indexing will not work"<<endl;
		    
		}
		chriRead = arr->chri;
		if(!gw->writeAlleleRecord(arr)){
		    cerr<<"GlacCAT: error writing record "<<arr<<endl;
		    exit(1);
		}
	    }

	    if(unsafemode){
		char * buffer;
		while(true){
		    buffer = gp.fillBuffer(sizeRecords);
		    if(buffer == 0){
			break;
		    }
		    if(!gw->writeBuffer(buffer,sizeRecords)){
			cerr<<"GlacCAT: error writing record "<<endl;
			exit(1);
		    }
		    delete(buffer);
		}
	    }else{
		while(gp.hasData()){
		    arr = gp.getData();
		    if(!gw->writeAlleleRecord(arr)){
			cerr<<"GlacCAT: error writing record "<<arr<<endl;
			exit(1);
		    }
		}
	    }
	}
    }

    delete(gw);

    return 0;
}

