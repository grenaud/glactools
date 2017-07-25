/*
 * vcf2glf
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "vcf2glf.h"

#define MIN(a,b) (((a)<(b))?(a):(b))


vcf2glf::vcf2glf(){

}

vcf2glf::~vcf2glf(){

}

string vcf2glf::usage() const{

    
    return string(string("") +"vcf2glf <options> [vcf file] [name sample] "+"\n"+
		  "\nThis program convert VCF files into glf (prints to the stdout)\n"+                       
		  //"\t"+"--bytes [#]" +"\t\t"+"Use either 2 or 3 bytes for storing allele count  (default: "+stringify(bytesForAC)+")\n"+
		  "\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
		  "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
		  
		  "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
		  "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
		  "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
		  "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
		  "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
		  
		  // "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
		  // "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 
		  
		  "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
		  "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
		  "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
		  "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
		  "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n");
    
    //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+


}


int vcf2glf::run(int argc, char *argv[]){
    //cout<<"ST "<<string(argv[0])<<endl;
    //last arg is program name
    for(int i=1;i<(argc-3);i++){ 

        if(string(argv[i]) == "-u"){
	    uncompressed=true;
            continue;
        }

        if(string(argv[i]) == "--fai"){
            fastaIndex=string(argv[i+1]);
            i++;
            continue;
        }

        // if( string(argv[i]) == "--bytes"  ){
        //     bytesForAC=destringify<int>(argv[i+1]);
        //     i++;
        //     continue;
	// }

        // if( string(argv[i]) == "--minPL"  ){
        //     minPLdiffind=destringify<int>(argv[i+1]);
	//     //            specifiedPL  =true;
        //     i++;
        //     continue;
	// }

        if(string(argv[i]) == "--minMap"){
            minMapabilitycutoff=destringify<double>(argv[i+1]);
            i++;
            continue;
	}
                 
	if(string(argv[i]) == "--minCov"){
	    minCovcutoff=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

       
	if(string(argv[i]) == "--maxCov"){
	    maxCovcutoff=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--minGQ"){
	    minGQcutoff=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

	if(string(argv[i]) == "--minMQ"){
	    minMQcutoff=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }
	

	if(string(argv[i]) == "--allowindel"){
	    allowCloseIndelProx =true;
	    continue;
	}

	if(string(argv[i]) == "--allowrm"){
	    allowRepeatMasking   =true;
	    continue;
	}

	if(string(argv[i]) == "--allowSysErr"){
	    allowSysErr     =true;
	    continue;
	}

	if(string(argv[i]) == "--allowall"){
	    allowall       =true;
	    continue;
	}

	if(string(argv[i]) == "--allowallMQ"){
	    allowallMQ      =true;
	    continue;
	}

	// if(strcmp(argv[i]) == "--useanc"){
	//     ancAllele       =true;
	//     continue;
	// }

	cerr<<"Wrong option "<<argv[i]<<endl;
	exit(1);
    }

    if(fastaIndex.size()==0){
	cerr<<"Must specify fai file "<<endl;
	exit(1);	
    }

    filtersVCF= new SetVCFFilters (minGQcutoff          ,
				   minMQcutoff          ,
				   minMapabilitycutoff  ,
				   !allowCloseIndelProx ,
				   !allowRepeatMasking  ,
				   !allowSysErr         ,
				   minCovcutoff      ,
				   maxCovcutoff  ,
				   allowall,
				   allowallMQ);

    cerr<<"Filter used: "<<*filtersVCF<<endl;
    string tempvcf_ = string(argv[argc-2]);
    //cerr<<"VCF file "<<tempvcf_<<endl;
    VCFreader vcfr   (tempvcf_,5);
    string namePop  = string(argv[argc-1]);
    // string epoFile  = string(argv[argc-1]);
    // string epoFileidx = epoFile+".tbi";
    // if(bytesForAC != 2 && bytesForAC != 3){
    // 	cerr<<"Specify either 2 or 3 bytes for the allele count"<<endl;
    // }

    cerr<<"Name pop "<<namePop<<endl;
    //cerr<<"EPO file "<<(string(argv[argc-1]))<<endl;

    // string epoFileidx  = epoFileidx+".tbi";
    // if(!isFile(epoFileidx)){
    // 	cerr<<"Tabix file epoFileidx not found"<<endl;
    // 	return 1;
    // }
    // igzstream epoFileFP;
    // epoFileFP.open(epoFile.c_str(), ios::in);    // open the streams
    // string epoLine;
    // string epoChr;
    // unsigned int epoCoord;
    //char epoREF;

    // if (epoFileFP.good()) {
    // 	//fine
    // }else{
    // 	cerr<<"Unable to open the file "<<epoFile<<endl;
    // 	exit(1);
    // }


    // getline(epoFileFP,epoLine);
    //read first line
    // if(1){
    // 	vector<string> fieldsEPO= allTokens(epoLine,'\t');
    // 	epoChr=fieldsEPO[0];
    // 	epoCoord=string2uint(fieldsEPO[1]);	
    // }


    // ReadTabix * rtEPO ;
    // string lineFromEPO;
    // bool lineLeftEPO;
    // bool cpgEPO=false;
    // bool firstLine=true;
    // char allel_chimp;
    // char allel_anc;

    // 75060,
    // 75070,
    // 5);
    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    readFastaIndex(fastaIndex,chrFound,genomeLength);


    
    string bgzf_file = "/dev/stdout";
    BGZF * fpBGZF;
    if(!uncompressed){
	fpBGZF = bgzf_open(bgzf_file.c_str(), "w");
    

	
    if (fpBGZF == NULL) { // region invalid or reference name not found
    	cerr<<"Cannot write to "<<bgzf_file<<endl;
    	exit(1);
    }else{
    	
    }
    }


    char bammagicstr [4] = {'B','A','M','\1'};    
    //magicstr="MST";
    if(uncompressed){
	if(write(1,&bammagicstr,sizeof(bammagicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
    }else{
	if( bgzf_write(fpBGZF, &bammagicstr,sizeof(bammagicstr)) != sizeof(bammagicstr) ){   cerr<<"Write error"<<endl;            return 1;   }     
    }
    char magicstr [5] = {'G','L','F', char(bytesForGL) ,'\1'};
    
    //magicstr="MST";
    //if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
    if(uncompressed){
	if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   }     
    }else{
	if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;            return 1;   }     
    }


    //HEADER

    string header="";
    //cout<<"#GLF"<<endl;    
    header+="#GLF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    //cout<<"#PG:"<<programLine<<endl;
    header+="#PG:"+programLine+"\n";
    //cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    header+="#GITVERSION: "+returnGitHubVersion(argv[0],"")+"\n";
    //cout<<"#DATE: "<<getDateString()<<endl;
    header+="#DATE: "+getDateString()+"\n";


    

    ///cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<namePop<<endl;

    //cout<<"#chr\tcoord\tREF,ALT\troot\t"<<namePop<<endl;


    map<string,uint16_t> chr2index;
    uint16_t     chrCurrentIndex=0;

    for(unsigned j=0;j<(chrFound.size());j++){
        //cout<<"#SQ\t"<<"SN:"<<chrFound[j].name<<"\tLN:"<<chrFound[j].length<<endl;
        //toflush<<"#SQ\t"<<"SN:"<<chrFound[j].name<<"\tLN:"<<chrFound[j].length<<endl;
        //toflush+=mp.getHeaderNoDefline();
	//header+= "#SQ\t"+"SN:"+chrFound[j].name+"\tLN:"+stringify(chrFound[j].length)+"\n";
	header+= string("#SQ\tSN:")+chrFound[j].name+"\tLN:"+stringify(chrFound[j].length)+"\n";
        chr2index[chrFound[j].name]=chrCurrentIndex++;
        if(chrCurrentIndex == 0xFFFF){
            cerr<<"Too many chromosomes for this build, more than 65535"<<endl;
            return 1;
        }
    }
    header+="#chr\tcoord\tREF,ALT\troot\tanc\t"+namePop+"\n";

    //toflush<<mp.getDefline()<<endl;
    
    uint32_t sizeHeader= header.size();
   
    if(uncompressed){
	if(write(1,&sizeHeader,sizeof(sizeHeader)) == -1 ){   cerr<<"Write error"<<endl;        return 1;   } 
    }else{
	if( bgzf_write(fpBGZF, &sizeHeader,sizeof(sizeHeader)) != sizeof(sizeHeader) ){   cerr<<"Write error"<<endl;            return 1;   }     
    }

    //char  header_cptr [sizeHeader];

    
    //strcpy(header_cptr,header.str().c_str());
    //cerrr<<
    for(uint32_t i=0;i<sizeHeader;i++){
        //char towrite=char(header_cptr[i]);
	char towrite=char(header[i]);
	if(uncompressed){
	    if(write(1,&towrite,sizeof(towrite)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
	}else{
	    if( bgzf_write(fpBGZF, &towrite,sizeof(towrite)) != sizeof(towrite) ){   cerr<<"Write error"<<endl;            return 1;   }     
	}
    }
    uint32_t sizePops=1;

    //if(write(1,&sizePops,sizeof(sizePops)) == -1 ){   cerr<<"Write error"<<endl;            return 1;   } 
    if(uncompressed){
	if(write(1,&sizePops,sizeof(sizePops)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
    }else{
	if( bgzf_write(fpBGZF, &sizePops,sizeof(sizePops)) != sizeof(sizePops) ){   cerr<<"Write error"<<endl;            return 1;   }     
    }

    


    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();


	
	if(passedFilters(toprint,filtersVCF)){



	    

	    
	    //pair<int,int> pairCount= toprint->returnLikelyAlleleCountForRefAlt(minPLdiffind);
	    //cerr<<pairCount.first<<"\t"<<pairCount.second<<endl;
	    //if(pairCount.first != 0 || pairCount.second != 0 ){
		char alt=(toprint->getAlt()=="."?'N':toprint->getAlt()[0]);
		string chimpString;
		string ancString;
		
		// cout<<allel_chimp<<"\t"<<toprint->getRef()[0]<<endl;
		

		//unresolved ancestral allele (A,C,G,T)
		//if(!isResolvedDNA(allel_chimp)){
		chimpString="0,0:0";					


		//if(!isResolvedDNA(allel_anc)){
		ancString="0,0:0";					

		
		if(chr2index.find(toprint->getChr()) == chr2index.end()){
		    cerr<<"Cannot find chr "<<toprint->getChr()<<" in index "<<endl;
		    return 1;
		}
		
		//write index of chr 2 bytes
		
		//if(write(1,&chr2index[toprint->getChr()],sizeof(chr2index[toprint->getChr()])) == -1 ) {   cerr<<"Write error"<<endl;         return 1;   }
		if(uncompressed){
		    if(write(1,&chr2index[toprint->getChr()],sizeof(chr2index[toprint->getChr()])) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &chr2index[toprint->getChr()],sizeof(chr2index[toprint->getChr()])) != sizeof(chr2index[toprint->getChr()]) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}

		
		//write coordinate 4 bytes
		uint32_t c = uint32_t(toprint->getPosition());
		//if(write(1,&c,    sizeof(c))  == -1 )    {   cerr<<"Write error"<<endl;         return 1;   }
		if(uncompressed){
		    if(write(1,&c,sizeof(c)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &c,sizeof(c)) != sizeof(c) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}
		

		uint8_t  tempCh;
		//uint16_t tempSh;

		//ref 1 byte
		tempCh= uint8_t(base2int(toprint->getRef()[0]));
		//if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return 1;   }
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}

		//alt 1 byte
		tempCh= uint8_t(base2int(alt));
		//if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return 1;   }
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}



		for(int i=0;i<2;i++){
		    //1 byte RR
		    tempCh= uint8_t( 0 );
		    //if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return 1;   }
		    if(uncompressed){
			if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		    }else{
			if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		    }


		    //1 byte RA
		    tempCh= uint8_t( 0 );
		    //if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return 1;   }
		    if(uncompressed){
			if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		    }else{
			if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		    }

		    //1 byte
		    tempCh= uint8_t( 0 );
		    //if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		    if(uncompressed){
			if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		    }else{
			if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		    }

		    //1 byte
		    tempCh= uint8_t( 0 );
		    //if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		    if(uncompressed){
			if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		    }else{
			if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		    }


		}
		
		//1 byte RR
		tempCh= uint8_t( MIN(toprint-> getPLHomoRef(),255) );
		//if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return 1;   }
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}


		//1 byte RA
		tempCh= uint8_t( MIN(toprint->getPLHetero(),255) );
		//if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return 1;   }
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}

		//1 byte
		tempCh= uint8_t( MIN(toprint->getPLHomoAlt(),255) );
		//if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}

		//1 byte
		tempCh= uint8_t( toprint->isCpg() );
		//if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return 1;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return 1;   }  
		}

		
		// for(unsigned j=0;j<(mp.getPopulationsNames()->size());j++){
		//     //2 bytes
		//     tempSh= short( dataRow->vectorAlleles->at(j).getRefCount() );
		//     if(write(1,&tempSh,sizeof(tempSh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		//     //2 bytes
		//     tempSh= short( dataRow->vectorAlleles->at(j).getAltCount() );
		//     if(write(1,&tempSh,sizeof(tempSh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		//     //1 byte
		//     tempCh= short( dataRow->vectorAlleles->at(j).getIsCpg() );
		//     if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;       return 1;   }
		// }
		
		// cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
		//     toprint->getRef()<<","<<
		//     alt<<"\t"<<
		//     chimpString<<"\t"<<
		//     ancString<<"\t"<<
		//     pairCount.first<<","<<pairCount.second<<":"<<(toprint->isCpg()?"1":"0")<<endl;	
		// }
	    //<<endl;
	}

    }

    if(!uncompressed){
	if(bgzf_close(fpBGZF) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;            return 1;   }  
    }

    cerr<<rejectFiltersTally()<<endl;
    //cerr<<*filtersVCF<<endl;
    //epoFileFP.close();
    delete(filtersVCF);
    //delete(rtEPO);
    cerr<<"Program terminated gracefully"<<endl;
    return 0;

}
