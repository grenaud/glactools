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
        
    return string(string("") +"glactools view [options] <gl|ac file> [region] "+
                  "\nThis program can view acf/glf files (prints to the stdout)\n"+                       
		  "\n"+
		  "\t"+"-b" + "\t\t\t"+"Produce binary compressed output (default: "+booleanAsString(uncompressed)+")\n"+
                  "\t"+"-u" + "\t\t\t"+"Produce binary but uncompressed output (default: "+booleanAsString(uncompressed)+")\n"+
		  "\n"+
                  "\t"+"" + ""+"For text output:\n"+		  
                  "\t"+"-h" + "\t\t\t"+"Produce defline     (default: "+booleanAsString(printdefline)+")\n"+
                  "\t"+"-H" + "\t\t\t"+"Produce full header (default: "+booleanAsString(printheader)+")\n");

}


int GlacViewer::run(int argc, char *argv[]){
    int lastOpt=1;

    //for(int i=1;i<(argc);i++){ 

    

    for(int i=1;i<(argc);i++){ 
	//cout<<i<<"\t"<<string(argv[i])<<endl;
	if((string(argv[i]) == "-") &&
	   ( i==(argc-1) || i==(argc-2) )
	   ){
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
	
	if(string(argv[i]) == "-b"){
	    printBin=true;
            continue;
        }
	
	if(string(argv[i]) == "-h"){
	    printdefline=true;
            continue;
        }

	if(string(argv[i]) == "-H"){
	    printheader=true;
            continue;
        }

        cerr<<"Error: unknown option "<<string(argv[i])<<endl;
        return 1;	
    }
    
    if(printBin && uncompressed){
        cerr<<"Error: cannot use both -u and -b"<<endl;
        return 1;	
    }

    if(printBin || uncompressed){
	if(printheader || printdefline){
	    cerr<<"Warning: if you specify either -u and -b, it will print the header automatically"<<endl;
	    //return 1;		    
	}
    }
    
    if(lastOpt == (argc-1)){//no region given

	string glacfile  = string(argv[lastOpt]);
	//cerr<<"glacfile "<<glacfile<<endl;
	//todo region retrieve
	GlacParser gp (glacfile);
	AlleleRecords * ar;
		

	GlacWriter * gw=NULL;

	if(uncompressed || printBin){//if binary
	    gw = new GlacWriter(gp.getSizePops(),
				gp.isGLFormat(),
				gp.isACFormat()?2:1,
				uncompressed);
	    if(!gw->writeHeader(gp.getHeader())){
		cerr<<"GlacViewer: error writing header "<<endl;
		exit(1);
	    }

	}else{//end if binary	
	    if(printheader)
		cout<<gp.getHeader()<<endl;
	    if(printdefline)
		cout<<gp.getDefline()<<endl;
	}

	
	while(gp.hasData()){

	    ar = gp.getData();

	    if(uncompressed || printBin){//if binary
		if(!gw->writeAlleleRecord(ar)){
		    cerr<<"GlacViewer: error record "<<*ar<<endl;
		    exit(1);
		}

	    }else{
		cout<<*ar<<endl;
	    }	    
	}

	delete(gw);

	return 0;
    }else{//else region given

	if(lastOpt == (argc-2)){//no region given
	    string glacfile  = string(argv[lastOpt]);
	    string region    = string(argv[lastOpt+1]);
	    //cout<<region<<endl;
	    unsigned int indexColon=0;
	    unsigned int indexDash =0;		
	    string chrName;
	    int start=-1;
	    int end=-1;
	    bool justChr;
	    for(unsigned int k=0;k<region.size();k++){
		if(region[k]=='-' && indexColon!=0)
		    indexDash=k;
		if(region[k]==':')
		    indexColon=k;
	    }
	    //cout<<indexColon<<endl;
	    if(indexColon == 0  ){
		//just chr
		justChr=true;
		chrName=region;
		
	    }else{
		justChr=false;
		chrName = region.substr(0,indexColon);
		start   = destringify<int>( region.substr(indexColon+1,indexDash-indexColon-1));
		end     = destringify<int>( region.substr(indexDash+1));
	    }
	    end = end+1;//I don't know why but this seems to work to give inclusive bounds
	    //cout<<justChr<<endl;
	    GlacParser gp (glacfile,glacfile+".bai",chrName, start, end, justChr);


	    GlacWriter * gw=NULL;

	    if(printBin || uncompressed){
		gw = new GlacWriter(gp.getSizePops(),
				    gp.isGLFormat(),
				    gp.isACFormat()?2:1,
				    uncompressed);
		if(!gw->writeHeader(gp.getHeader())){
		    cerr<<"GlacViewer: error writing header "<<endl;
		    exit(1);
		}
	    }else{
		if(printheader)
		    cout<<gp.getHeader()<<endl;
		if(printdefline)
		    cout<<gp.getDefline()<<endl;
	    }


	    AlleleRecords * ar;
	    while(gp.hasData()){

		ar = gp.getData();
		if(uncompressed || printBin){//if binary
		    if(!gw->writeAlleleRecord(ar)){
			cerr<<"GlacViewer: error record "<<*ar<<endl;
			exit(1);
		    }

		}else{
		    cout<<*ar<<endl;
		}	    
		
	    }

	    delete(gw);

	    return 0;
	}else{
	    cerr<<"GlacViewer: specify file [reg] or just the file"<<endl;
	    return 1;
	}
	return 1;
    }
    return 0;
}

