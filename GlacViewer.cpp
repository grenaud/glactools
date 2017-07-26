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

                  "\t"+"-h" + "\t\t\t"+"Produce defline     (default: "+booleanAsString(printdefline)+")\n"+
                  "\t"+"-H" + "\t\t\t"+"Produce full header (default: "+booleanAsString(printheader)+")\n");

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
	//todo region retrieve
	GlacParser gp (glacfile);
	AlleleRecords * test;

	string bgzf_file = "/dev/stdout";
	BGZF * fpBGZF    = NULL;
	if(printBin && !uncompressed){
	    fpBGZF = bgzf_open(bgzf_file.c_str(), "w");
	    	    	
	    if (fpBGZF == NULL) { // region invalid or reference name not found
		cerr<<"Cannot write to "<<bgzf_file<<endl;
		exit(1);
	    }else{		
	    }
	}

	size_t sizeRecord ;

	if(uncompressed || printBin){//if binary
	    sizeRecord = gp.getSizeRecord();
	    char bammagicstr [4] = {'B','A','M','\1'};    

	    if(uncompressed){
		if(write(1,&bammagicstr,sizeof(bammagicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
	    }else{
		if( bgzf_write(fpBGZF, &bammagicstr,sizeof(bammagicstr)) != sizeof(bammagicstr) ){   cerr<<"Write error"<<endl;   return 1;   }     
	    }

	    if(gp.isACFormat()){
		char magicstr [5] = {'A','C','F', char(bytesForAC) ,'\1'};
		if(uncompressed){
		    if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
		}else{
		    if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;  return 1;   }     
		}
	    }else{		
		char magicstr [5] = {'G','L','F', char(bytesForGL) ,'\1'};   
		if(uncompressed){
		    if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
		}else{
		    if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;  return 1;   } 
		}
	    }
	    
	    uint32_t sizeHeader= gp.getHeader().size();
	    if(uncompressed){
		if(write(1,&sizeHeader,sizeof(sizeHeader)) == -1 ){   cerr<<"Write error"<<endl;        return 1;   } 
	    }else{
		if( bgzf_write(fpBGZF, &sizeHeader,sizeof(sizeHeader)) != sizeof(sizeHeader) ){   cerr<<"Write error"<<endl; return 1;   }     
	    }

	    for(uint32_t i=0;i<sizeHeader;i++){
		char towrite=char(gp.getHeader()[i]);
		if(uncompressed){
		    if(write(1,&towrite,sizeof(towrite)) == -1 ){   cerr<<"Write error"<<endl;  return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &towrite,sizeof(towrite)) != sizeof(towrite) ){   cerr<<"Write error"<<endl; return 1;   }     
		}
	    }

	    uint32_t sizePops=gp.getSizePops();
	    if(uncompressed){
		if(write(1,&sizePops,sizeof(sizePops)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
	    }else{
		if( bgzf_write(fpBGZF, &sizePops,sizeof(sizePops)) != sizeof(sizePops) ){   cerr<<"Write error"<<endl;            return 1;   }     
	    }
	}else{//end if binary	
	    if(printheader)
		cout<<gp.getHeader()<<endl;
	    if(printdefline)
		cout<<gp.getDefline()<<endl;
	}
	//TODO stopped here wednesday
	char buffer [ gp.getSizeRecord() ];

	
	while(gp.hasData()){

	    test = gp.getData();
	    // cout<<"run"<<endl;
	    // cout<<test<<endl;
	    cout<<*test<<endl;
	    
	}
	
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
	    //cout<<justChr<<endl;
	    GlacParser gp (glacfile,glacfile+".bai",chrName, start, end, justChr);

	    string bgzf_file = "/dev/stdout";
	    BGZF * fpBGZF    = NULL;
	    if(printBin && !uncompressed){
		fpBGZF = bgzf_open(bgzf_file.c_str(), "w");
	    	    	
		if (fpBGZF == NULL) { // region invalid or reference name not found
		    cerr<<"Cannot write to "<<bgzf_file<<endl;
		    exit(1);
		}else{		
		}
	    }

	    size_t sizeRecord ;

	    if(uncompressed || printBin){//if binary
		sizeRecord = gp.getSizeRecord();
		char bammagicstr [4] = {'B','A','M','\1'};    

		if(uncompressed){
		    if(write(1,&bammagicstr,sizeof(bammagicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
		}else{
		    if( bgzf_write(fpBGZF, &bammagicstr,sizeof(bammagicstr)) != sizeof(bammagicstr) ){   cerr<<"Write error"<<endl;   return 1;   }     
		}

		if(gp.isACFormat()){
		    char magicstr [5] = {'A','C','F', char(bytesForAC) ,'\1'};
		    if(uncompressed){
			if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
		    }else{
			if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;  return 1;   }     
		    }
		}else{		
		    char magicstr [5] = {'G','L','F', char(bytesForGL) ,'\1'};   
		    if(uncompressed){
			if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
		    }else{
			if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;  return 1;   } 
		    }
		}
	    
		uint32_t sizeHeader= gp.getHeader().size();
		if(uncompressed){
		    if(write(1,&sizeHeader,sizeof(sizeHeader)) == -1 ){   cerr<<"Write error"<<endl;        return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &sizeHeader,sizeof(sizeHeader)) != sizeof(sizeHeader) ){   cerr<<"Write error"<<endl; return 1;   }     
		}

		for(uint32_t i=0;i<sizeHeader;i++){
		    char towrite=char(gp.getHeader()[i]);
		    if(uncompressed){
			if(write(1,&towrite,sizeof(towrite)) == -1 ){   cerr<<"Write error"<<endl;  return 1;   } 
		    }else{
			if( bgzf_write(fpBGZF, &towrite,sizeof(towrite)) != sizeof(towrite) ){   cerr<<"Write error"<<endl; return 1;   }     
		    }
		}

		uint32_t sizePops=gp.getSizePops();
		if(uncompressed){
		    if(write(1,&sizePops,sizeof(sizePops)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &sizePops,sizeof(sizePops)) != sizeof(sizePops) ){   cerr<<"Write error"<<endl;            return 1;   }     
		}
	    }else{//end if binary	

		if(printheader)
		    cout<<gp.getHeader()<<endl;
		if(printdefline)
		    cout<<gp.getDefline()<<endl;
	    }

	    AlleleRecords * test;
	    while(gp.hasData()){

		test = gp.getData();
		// cout<<"run"<<endl;
		// cout<<test<<endl;
		cout<<*test<<endl;
		//cout<<"GlacViewer return= "<<*test<<endl;
	    
	    }

	    //cout<<"l"<<endl;
	    return 0;
	}else{
	    cerr<<"GlacViewer: specify file [reg] or just the file"<<endl;
	    return 1;
	}
	return 1;
    }
    return 0;
}

