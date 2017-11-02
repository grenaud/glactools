/*
 * Vcf2ACF
 * Date: Jul-23-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "Vcf2ACF.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

//#define DEBUGPOS 9484430

Vcf2ACF::Vcf2ACF(){

}

Vcf2ACF::~Vcf2ACF(){

}


void Vcf2ACF::setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO){

    lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
    if(!lineLeftEPO){
	cerr<<"Error, missing data in the EPO file"<<endl;
	exit(1);
    }

    vector<string> fieldsEPO  = allTokens(lineFromEPO,'\t');
    epoChr                   = fieldsEPO[0];
    epoCoord                 = string2uint(fieldsEPO[1]);					
    if(fieldsEPO[9] == "1")
	cpgEPO=true;		    
    else
	cpgEPO=false;		    



    allel_anc   = fieldsEPO[3][0];//inferred ancestor
    allel_chimp = fieldsEPO[4][0];//chimp;


}

string Vcf2ACF::usage() const{

    
    return string(string("") +"vcf2acf [options] <vcf file> <name sample> "+"\n"+
		  "\nThis program converts a VCF file with a single sample into acf (prints to the stdout)\n"+                       
		  // "\nThe name of the samples has to be comma-separated\n"+                       
		  // "\ne.g. vcf2acf myfile.vcf.gz Ind1,Ind2,Ind10\n"+                       
		  
		  //"\t"+"--bytes [#]" +"\t\t"+"Use either 2 or 3 bytes for storing allele count  (default: "+stringify(bytesForAC)+")\n"+
		  "\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
		  "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+

		  "\t"+"--epo [EPO file]"       +"\t" +"Use file as EPO alignment to set the (default: none)\n"+   
		  "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   
		  
		  "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
		  "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
		  "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
		  "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
		  "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
		  "\t"+"--onlyGT"        +"\t\t" +"Do not use PL values for alleles, simply use genotypes (GT)      (default: "+booleanAsString(onlyGT)+")\n"+ 
		  "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
		  // "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 
		  
		  "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
		  "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
		  "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
		  "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
		  "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n");
    
    //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+


}


int Vcf2ACF::run(int argc, char *argv[]){

    int lastOpt=1;
    bool specifiedPL  = false;
    
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


	if(strcmp(argv[i],"--epo") == 0 ){
            epoFile=string(argv[i+1]);
	    epoFileB=true;
            i++;
            continue;
	}

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

        if( string(argv[i]) == "--onlyGT"  ){
	    onlyGT  = true;
	    //            specifiedPL  =true;
            continue;
	}

        if( string(argv[i]) == "--minPL"  ){
            minPLdiffind = destringify<int>(argv[i+1]);
	    specifiedPL  = true;
            i++;
            continue;
	}

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

	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    if(specifiedPL  && onlyGT){
	cerr<<"Cannot both operate on GT and PL simultaneously "<<endl;
	return 1;	
    }
    
    if(fastaIndex.size()==0){
	cerr<<"Must specify fai file "<<endl;
	return 1;	
    }

    if(lastOpt != (argc-2)){
	cerr<<"The last 2 arguments are <vcf file> <name sample> "<<endl;
	return 1;		
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
    string tempvcf_ = string(argv[lastOpt]);
    //cerr<<"VCF file "<<tempvcf_<<endl;
    VCFreader vcfr   (tempvcf_,5);
    string namePop  = string(argv[lastOpt+1]);
    //cerr<<"namePop "<<namePop<<endl;
    //vector<string> namePopToReturn = allTokens(namePop,',');
    // string epoFile  = string(argv[argc-1]);
    string epoFileidx = epoFile+".tbi";
    if(bytesForAC != 2 && bytesForAC != 3){
	cerr<<"Specify either 2 or 3 bytes for the allele count"<<endl;
    }

    //cerr<<"Name pop "<<namePop<<endl;
    //cerr<<"EPO file "<<(string(argv[argc-1]))<<endl;

    // string epoFileidx  = epoFileidx+".tbi";
    // if(!isFile(epoFileidx)){
    // 	cerr<<"Tabix file epoFileidx not found"<<endl;
    // 	return 1;
    // }
    // igzstream epoFileFP;
    // epoFileFP.open(epoFile.c_str(), ios::in);    // open the streams
    // string epoLine;
    string epoChr;
    unsigned int epoCoord;
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


    ReadTabix * rtEPO =NULL;
    string lineFromEPO;
    bool lineLeftEPO;
    bool cpgEPO=false;
    bool firstLine=true;
    char allel_chimp;
    char allel_anc;

    // 75060,
    // 75070,
    // 5);
    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    readFastaIndex(fastaIndex,chrFound,genomeLength);

    GlacWriter * gw=NULL;

    gw = new GlacWriter(1,    //gp.getSizePops(),
			false, //gp.isGLFormat(),
			2,//gp.isACFormat()?2:1,
			1,//compression threads
			uncompressed);
    

	




    //HEADER

    string header="";
    //cout<<"#ACF"<<endl;    
    header+="#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    //cout<<"#PG:"<<programLine<<endl;
    header+="#PG:"+programLine+"\n";
    header+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";
    //cout<<"#DATE: "<<getDateString()<<endl;
    header+="#DATE: "+getDateString()+"\n";
    header+="#VCF2ACF:\n";


    


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



    
    if(!gw->writeHeader(header)){
	cerr<<"GlacViewer: error writing header "<<endl;
	return 1;
    }
    


    


    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();

#ifdef DEBUGPOS
	bool debugPosition=false;
	if(toprint->getPosition() == DEBUGPOS){
	    debugPosition=true;
	    cerr<<*toprint<<endl;
	}
#endif
	if(passedFilters(toprint,filtersVCF)){


#ifdef DEBUGPOS
	    if(debugPosition)
		cerr<<"passed "<<endl;
#endif
	    
	    if(firstLine){
		firstLine=false;
		if(epoFileB){
		    rtEPO = new ReadTabix( epoFile.c_str()  , 
					   epoFileidx.c_str()  , 
					   toprint->getChr(), 
					   int(toprint->getPosition()),INT_MAX ); //the destructor should be called automatically
		    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		}
		// lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
		// vector<string> fieldsEPO = allTokens(lineFromEPO,'\t');
		// epoChr                   = fieldsEPO[0];
		// epoCoord                 = string2uint(fieldsEPO[1]);	
                // if(fieldsEPO[9] == "1")
		//     cpgEPO=true;
		// else
		//     cpgEPO=false;
		// if(ancAllele){
		//     allel_chimp = fieldsEPO[3][0];//inferred ancestor
		// }else{
		//     allel_chimp = fieldsEPO[4][0];//chimp;
		// }		    
	    }

	    if(epoFileB){

		if(!lineLeftEPO){
		    cerr<<"Error, no data in the EPO file "<< toprint->getChr() <<":"<< int(toprint->getPosition()) <<endl;
		    return 1;
		}
	    
		if(epoChr != toprint->getChr()){
		    //cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<toprint->getChr()<<endl;
		    //return 1;

		    rtEPO->repositionIterator(toprint->getChr()  , int(toprint->getPosition()) ,INT_MAX);
		    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		    if(epoChr != toprint->getChr()){
			cerr<<"Error, the repositioning did not work, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<toprint->getChr()<<endl;
			return 1;
		    }

		}
	    
	    

		while(epoCoord != toprint->getPosition()){
		    if(epoCoord > toprint->getPosition()){
			cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(*toprint)<<"\t"<<lineFromEPO<<endl;
			return 1;
		    }

		    if( (toprint->getPosition() - epoCoord ) >= limitToReOpenFP){ //seeking in the file
			rtEPO->repositionIterator(toprint->getChr() , int(toprint->getPosition()),INT_MAX);
		    }


		    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		    // lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
		    // vector<string> fieldsEPO = allTokens(lineFromEPO,'\t');
		    // epoChr                   = fieldsEPO[0];
		    // epoCoord                 = string2uint(fieldsEPO[1]);					
		    // if(fieldsEPO[9] == "1")
		    //     cpgEPO=true;		    
		    // else
		    //     cpgEPO=false;		    
		    // if(ancAllele){
		    //     allel_chimp = fieldsEPO[3][0];//inferred ancestor
		    // }else{
		    //     allel_chimp = fieldsEPO[4][0];//chimp;
		    // }

		    // if(!lineLeftEPO){
		    //     cerr<<"Error, missing data in the EPO file"<<*toprint<<endl;
		    //     return 1;
		    // }
		}
	    
		if(epoCoord != toprint->getPosition()){
		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(*toprint)<<"\t"<<lineFromEPO<<endl;
		    return 1;
		}
	    }
	    

	    
#ifdef DEBUGPOS
	    if(debugPosition)
		cerr<<"debug pos"<<endl;
		
#endif

	    pair<int,int> pairCount;
	    if(onlyGT)
		pairCount = toprint->returnLikelyAlleleCountForRefAltJustGT();
	    else
		pairCount = toprint->returnLikelyAlleleCountForRefAlt(minPLdiffind);

#ifdef DEBUGPOS
	    if(debugPosition){
		cerr<<pairCount.first<<"\t"<<pairCount.second<<"\t"<<minPLdiffind<<endl;
		//return 1;
	    }
#endif
	    
	    
	    if(pairCount.first != 0 || pairCount.second != 0 ){
		char alt=(toprint->getAlt()=="."?'N':toprint->getAlt()[0]);
		// string chimpString;
		// string ancString;
		SingleAllele root;
		SingleAllele anc;

		// cout<<allel_chimp<<"\t"<<toprint->getRef()[0]<<endl;
		if(!epoFileB){ //no epo file
		    // chimpString="0,0:0";
		    // ancString="0,0:0";
		    root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
		    anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

		}else{	//epoFileB	
		    if(!isResolvedDNA(allel_chimp)){
			//chimpString="0,0:0";					
			root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
		    }
		    //resolved ancestral allele
		    else{
			if(allel_chimp == toprint->getRef()[0]){//no diff between chimp and reference
			    //chimpString="1,0:"+string(cpgEPO?"1":"0");
			    root.setRefCount(1); root.setAltCount(0);  root.setIsCpg(cpgEPO); 
			}else{
			    if(alt == 'N'){//no alt defined, the chimp becomes the alt			    
				alt = allel_chimp;
				//chimpString="0,1:"+string(cpgEPO?"1":"0");
				root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 
			    }else{
				if(alt == allel_chimp){//alt is chimp
				    //chimpString="0,1:"+string(cpgEPO?"1":"0");
				    root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 
				}else{ //tri-allelic site, discard
				    continue;
				}
			    }
			}
		    }	

		    if(!isResolvedDNA(allel_anc)){
			//ancString="0,0:0";
			anc.setRefCount(0); anc.setAltCount(0);  anc.setIsCpg(false); 
		    }
		    //resolved ancestral allele
		    else{
			if(allel_anc == toprint->getRef()[0]){//no diff between ancestor and reference
			    //ancString="1,0:"+string(cpgEPO?"1":"0");
			    anc.setRefCount(1); anc.setAltCount(0);  anc.setIsCpg(cpgEPO); 
			}else{
			    if(alt == 'N'){//no alt defined, the ancestor becomes the alt			    
				alt = allel_anc;
				//ancString="0,1:"+string(cpgEPO?"1":"0");			    
				anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 
			    }else{
				if(alt == allel_anc){//alt is ancestor
				    //ancString="0,1:"+string(cpgEPO?"1":"0");
				    anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 
				}else{ //tri-allelic site, discard
				    continue;
				}
			    }
			}
		    }
		}//end if we have epo file

		//unresolved ancestral allele (A,C,G,T)
		//if(!isResolvedDNA(allel_chimp)){
		//chimpString="0,0:0";					


		//if(!isResolvedDNA(allel_anc)){
		//ancString="0,0:0";					

		
		if(chr2index.find(toprint->getChr()) == chr2index.end()){
		    cerr<<"Cannot find chr "<<toprint->getChr()<<" in index "<<endl;
		    return 1;
		}
		
		
		AlleleRecords arToWrite (false);
		arToWrite.chri          = chr2index[toprint->getChr()];
		arToWrite.coordinate    = toprint->getPosition();
		arToWrite.sizePops      = 1;
		arToWrite.ref           = toprint->getRef()[0];		
		arToWrite.alt           = alt;		
		
		arToWrite.vectorAlleles = new vector<SingleAllele>  ();
		// SingleAllele root;
		// SingleAllele anc;
		SingleAllele sample (pairCount.first, pairCount.second, toprint->isCpg());

		//takes care of the case where for low coverage
		if(!root.hasAlt()    &&
		   !anc.hasAlt()    && 
		   !sample.hasAlt()   ){
		    arToWrite.alt           = 'N';
		}
		   
		arToWrite.vectorAlleles->push_back(root);
		arToWrite.vectorAlleles->push_back(anc);
		arToWrite.vectorAlleles->push_back(sample);
		
		if(!gw->writeAlleleRecord(&arToWrite)){
		    cerr<<"Vcf2ACF: error writing header "<<endl;
		    return 1;
		}


	    }else{
		//the allele count is 0
	    }
	    //<<endl;
	}else{

#ifdef DEBUGPOS
	    if(debugPosition)
		cerr<<"did not passed "<<endl;
#endif
	    
	}

    }

    // if(!uncompressed){
    // 	if(bgzf_close(fpBGZF) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;            return 1;   }  
    // }


    cerr<<rejectFiltersTally()<<endl;
    //cerr<<*filtersVCF<<endl;
    //epoFileFP.close();
    delete(filtersVCF);
    delete(gw);
    if(epoFileB)
	delete(rtEPO);
    cerr<<"Program terminated gracefully"<<endl;
    return 0;

}
