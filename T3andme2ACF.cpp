
#include "T3andme2ACF.h"


// #include "BAMTableObj.h"
// #include "BAMTABLEreader.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;



void T3andme2ACF::setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO){

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


    allel_ref   = fieldsEPO[2][0];//reference allele
    allel_anc   = fieldsEPO[3][0];//inferred ancestor
    allel_chimp = fieldsEPO[4][0];//chimp;


}



T3andme2ACF::T3andme2ACF(){

}

T3andme2ACF::~T3andme2ACF(){

}


string T3andme2ACF::usage() const{
    const string usage=string("")+ " 23andme2acf [options] <23andme file> <name sample>"+"\n"+
			      "\nThis program convert 23andme files into ACF (prints to the stdout)\n"+
			      
			      "\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
			      "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+

			      "\t"+"--epo [EPO file]"       +"\t" +"Use file as EPO alignment to set the (default: none)\n"+   
			      "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   
		  
			      "";
    return usage;
}

int T3andme2ACF::run(int argc, char *argv[]){

    // cout<<INT_MAX<<endl;
    // bool ancAllele           = false;

    int lastOpt=1;
			      

    //last arg is program name
    for(int i=1;i<(argc);i++){ 
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

	// if(strcmp(argv[i],"--useanc") == 0 ){
	//     ancAllele       =true;
	//     continue;
	// }


	cerr<<"Wrong option "<<argv[i]<<endl;
	exit(1);
    }


    if(lastOpt != (argc-2)){
	cerr<<"The last 2 arguments are <23andMe file> <name sample> "<<endl;
	return 1;		
    }


    if(fastaIndex.size()==0){
	cerr<<"Must specify fai file "<<endl;
	return 1;	
    }

    if(!epoFileB){
	cerr<<"Must specify EPO file as this file format does not have the reference information "<<endl;
	return 1;	
    }

    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    readFastaIndex(fastaIndex,chrFound,genomeLength);

    //cerr<<"Filter used: "<<*filtersVCF<<endl;
    //VCFreader vcfr   (string(argv[argc-3]),5);
    // BAMTABLEreader btr  (),5);
    string inputFile23 = string(argv[lastOpt]);
    igzstream myfile;
    myfile.open(inputFile23.c_str(), ios::in);

    if (!myfile.good()){
        cerr << "Unable to open file "<<inputFile23<<endl;
        return 1;
    }

    
    string namePop  = string(argv[lastOpt+1]);
    // string epoFile  = string(argv[argc-1]);
    string epoFileidx = epoFile+".tbi";
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
    char allel_ref;
    char allel_anc;

    // 75060,
    // 75070,
    // 5);
    //header
    unsigned int totalRec=0;
    unsigned int writtenRec=0;

    string header="";

    header+="#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    header+="#PG:"+programLine+"\n";;
    header+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";;

    header+="#DATE: "+getDateString()+"\n";;
    header+="#23ANDME2ACF:\n";
    map<string,uint16_t> chr2index;
    uint16_t     chrCurrentIndex=0;

    for(unsigned j=0;j<(chrFound.size());j++){
        header+= string("#SQ\tSN:")+chrFound[j].name+"\tLN:"+stringify(chrFound[j].length)+"\n";
        chr2index[chrFound[j].name]=chrCurrentIndex++;
        if(chrCurrentIndex == 0xFFFF){
            cerr<<"Too many chromosomes for this build, more than 65535"<<endl;
            return 1;
        }
    }



    header+="#chr\tcoord\tREF,ALT\troot\tanc\t"+namePop+"\n";
    GlacWriter * gw=NULL;

    gw = new GlacWriter(1,    //gp.getSizePops(),
                        false, //gp.isGLFormat(),
                        2,//gp.isACFormat()?2:1,
			1,//compression threads
                        uncompressed);
    

    if(!gw->writeHeader(header)){
	cerr<<"GlacViewer: error writing header "<<endl;
	exit(1);
    }
    



    string line;
    bool hasData= (bool)getline (myfile,line);
    while (hasData ){
	if(strBeginsWith(line,"#")){
	    hasData = (bool)getline (myfile,line);
	    continue;
	}
	 
	totalRec++;

	// while(btr.hasData()){
	// BAMTableObj * toprint=btr.getData();
	vector<string> token= allTokens(line,'\t');
	// cerr<<"line "<<line<<endl;
	string chr    = token[1];
	unsigned int pos = destringify<unsigned int>(token[2]);
	string genotype    = token[3];

	// if(passedFilters(toprint,filtersVCF)){

	if(firstLine){
	    firstLine=false;
	    if(epoFileB){

		rtEPO = new ReadTabix( epoFile.c_str()  , 
				       epoFileidx.c_str()  , 
				       chr,
				       int(pos),INT_MAX ); //the destructor should be called automatically

		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
	    }
	}

	if(epoFileB){

	    if(!lineLeftEPO){
		cerr<<"Error, no data in the EPO file "<< chr <<":"<< int(pos) <<endl;
		return 1;
	    }

	    if(epoChr != chr){
		//cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<chr<<endl;
		//return 1;
		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
		}	       	    
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
	    }


	    // cout<<"2"<<endl;
	    while(epoCoord != pos){
		if(epoCoord > pos){
		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(line)<<"\t"<<lineFromEPO<<endl;
		    return 1;
		}
		
		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
		}
		
	    
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);

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


	    if(epoCoord != pos){
		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<line<<"\t"<<lineFromEPO<<endl;
		return 1;
	    }
	}//end if epoFileB
	//int refIdx  =base2int(allel_ref);
	// int chimpIdx=base2int(allel_chimp);
	char alt='N';
	string s="ACGT";
	string chimpString;
	string ancString;
	SingleAllele root;
	SingleAllele anc;

	//need an allele that is A,C,G,T
	//cout<<pos<<"\t"<<allel_ref<<endl;
	if(!isResolvedDNA(allel_ref)){
	    goto nextline;
	}

	// cout<<"3"<<endl;

	//determine alternative allele
	for(int i=0;i<4;i++){
	    if(s[i] != allel_ref){
		if(genotype.find(s[i]) != string::npos){ //new non-ref allele found
		    //toprint->hasAllele(i+1) ){
		    alt=s[i];
		}
	    }
	}

	// cout<<*toprint<<endl;
	// cout<<lineFromEPO<<endl;
		
	// SingleAllele root;
	// SingleAllele anc;
	if(!epoFileB){ //no epo file
	    // chimpString="0,0:0";
	    // ancString="0,0:0";
	    root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
	    anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

	}else{
	    //unresolved ancestral allele (A,C,G,T)
	    if(!isResolvedDNA(allel_chimp)){
		//chimpString="0,0:0"; 
		root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
	    }
	    //resolved ancestral allele
	    else{
		if(allel_chimp == allel_ref){//no diff between chimp and reference
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
			    //continue;
			    goto nextline;
			}
		    }
		}
	    }


	    if(!isResolvedDNA(allel_anc)){
		//ancString="0,0:0"; 
		root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 
	    }
	    //resolved ancestral allele
	    else{
		if(allel_anc == allel_ref){//no diff between ancestor and reference
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
			    //continue;
			    goto nextline;
			}
		    }
		}
	    }
	}

	// cout<<"4"<<endl;



	if(alt!='N')
	    if(genotype.find(alt) == string::npos ){   // has alternative
		goto nextline;
	    }

	// cerr<<chr<<" "<<pos<<"alt "<<alt<<endl;
//if( (alt!='N' && genotype.find(alt) != string::npos) ){   // has alternative
	if(true){

	    int refCount=0;
	    int altCount=0;
	    string rr = stringify(allel_ref)+ stringify(allel_ref);
	    string ra = stringify(allel_ref)+ stringify(alt);
	    string ar = stringify(alt)      + stringify(allel_ref);
	    string aa = stringify(alt)      + stringify(alt);

	    if(genotype == rr ){
		refCount=2;
	    }else{
		if(genotype == aa){
		    altCount=2;
		}else{
		    if( (genotype == ra) ||
			(genotype == ar) ){
			refCount=1;
			altCount=1;
		    }else{
			//goto nextline;
			cerr<<"Error: Potential error in the 23 and me file where the reference allele is not there "<<line<<"\t"<<lineFromEPO<<endl; //error
			goto nextline;
		    }			
		}
	    }

	    //cout<<refIdx<<endl;
	    writtenRec++;
	    // cerr<<"WRITE "<<chr<<" "<<pos<<" alt "<<alt<<endl;


	    AlleleRecords arToWrite (false);
	    arToWrite.chri          = chr2index[chr];
	    arToWrite.coordinate    = pos;
	    arToWrite.sizePops      = 1;
	    arToWrite.ref           = allel_ref;
	    arToWrite.alt           = alt;		
		
	    arToWrite.vectorAlleles = new vector<SingleAllele>  ();
	    // SingleAllele root;
	    // SingleAllele anc;
	    SingleAllele sample (refCount, altCount, 0);//no CpG provided
		
	    arToWrite.vectorAlleles->push_back(root);
	    arToWrite.vectorAlleles->push_back(anc);
	    arToWrite.vectorAlleles->push_back(sample);
	    
	    if(!gw->writeAlleleRecord(&arToWrite)){
		cerr<<"T3ancme2ACF: error writing header "<<endl;
		exit(1);
	    }	    
	    // cout<<chr<<"\t"<< pos<<"\t"<<
	    // 	allel_ref<<","<<
	    // 	alt<<"\t"<<
	    // 	chimpString<<"\t"<<
	    // 	ancString<<"\t"<<
	    // 	refCount<<","<<
	    // 	altCount<<
	    // 	":"<<("0")<<endl;//no cpg info

	    // }else{
	    //     goto nextline;
	    //     //continue;//triallelic, skip
	    // }
	}
        nextline:
        
        hasData = (bool)getline (myfile,line);
	
    }

    //epoFileFP.close();
    //    delete(filtersVCF);
    if(rtEPO != NULL)
	delete(rtEPO);
    delete(gw);
    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRec<<" records, wrote "<<writtenRec<<" terminated gracefully"<<endl;

    return 0;
}

