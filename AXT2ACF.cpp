
#include "AXT2ACF.h"




using namespace std;



AXT2ACF::AXT2ACF(){

}

AXT2ACF::~AXT2ACF(){

}


string AXT2ACF::usage() const{
    const string usage=string("")+ " axt2acf [options]  <chr name> <name sample>  <axt file>\n"+

	"\nThis program will parse an AXT alignment for a single chromosome and output an ACF file (to STDOUT)\n"+
	"The first sequence has to be the genome of the reference. For instance, to import chimp for the human reference:\n"+
	"  glactools axt2acf  1 PanTro5 chr1.hg19.panTro5.net.axt.gz > output.acf.gz\n"
	"\n"+
	"\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+

		  
	"";
    return usage;
}

int AXT2ACF::run(int argc, char *argv[]){

    // cout<<INT_MAX<<endl;
    // bool ancAllele           = false;

    int lastOpt=1;
    string fastaIndex;

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

        if(string(argv[i]) == "-u"){
	    uncompressed=true;
            continue;
        }

        if(string(argv[i]) == "--fai"){
            fastaIndex=string(argv[i+1]);
            i++;
            continue;
        }


	cerr<<"Wrong option "<<argv[i]<<endl;
	exit(1);
    }


    if(lastOpt != (argc-3)){
	cerr<<"The last 3 arguments are <chr name> <name sample> <axt file>  "<<lastOpt<<" "<<argc<<endl;
	return 1;		
    }

    if(fastaIndex.size()==0){
	cerr<<"Must specify fai file "<<endl;
	return 1;	
    }

    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    readFastaIndex(fastaIndex,chrFound,genomeLength);

    //cerr<<"Filter used: "<<*filtersVCF<<endl;
    //VCFreader vcfr   (string(argv[argc-3]),5);
    // BAMTABLEreader btr  (),5);

    uint16_t chrIdx   = UINT16_MAX;
    string chrname    = string(argv[lastOpt+0]);
    string nameSample = string(argv[lastOpt+1]);
    string axtfile    = string(argv[lastOpt+2]);
    
    
    bool foundChr=false;
    for(unsigned int i=0;i<chrFound.size();i++){ 
	if(chrFound[i].name == chrname){
	    chrIdx   = uint16_t(i);
	    foundChr = true;
	    break;
	}
    }
    if(!foundChr){
        cerr << "Unable to find chromosome "<<chrname<<" in fasta index"<<endl;
        return 1;
    }

    
    unsigned int totalRec=0;
    unsigned int writtenRec=0;

    string header="";

    header+="#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    header+="#PG:"+programLine+"\n";
    header+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";

    header+="#DATE: "+getDateString()+"\n";
    header+="#AXT2ACF:"+chrname+" "+nameSample+" "+axtfile+"\n";

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



    header+="#chr\tcoord\tREF,ALT\troot\tanc\t"+nameSample+"\n";
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
    string lineSp1;
    string lineSp2;
    
    igzstream myAXTFile;
    myAXTFile.open(axtfile.c_str(), ios::in);

    SingleAllele root;
    SingleAllele anc;

    root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
    anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

    if (myAXTFile.good()){
	while ( getline (myAXTFile,line)){
	    if(strBeginsWith(line,"#"))
		continue;
	 
	    vector<string> allToks = allTokens(line,' ');

	    if(allToks[1]  != chrname){
		if(allToks[1]  != "chr"+chrname){		    
		    cerr<<"The chromosome name "<<allToks[1]<<" does not match the one provided on the command line ("<<"chr"<<chrname<<"), full line:"<<line<<endl;
		    return 1;
		}else{
		    chrname="chr"+chrname;
		}
	    }

	    unsigned int  startC = destringify<unsigned int>(allToks[2]);
	    unsigned int  endC   = destringify<unsigned int>(allToks[3]);
	    unsigned int  coordC =startC;
	    unsigned int  lastC = startC;
	    bool lastCharWasC=false;

	    if(allToks.size() == 9){
		getline (myAXTFile,lineSp1);
		getline (myAXTFile,lineSp2);

		//string  lastToPrintS="";
		AlleleRecords lastToPrint (false);
		bool          lastToPrintSet = false;

		for(unsigned i=0;i<lineSp1.size();i++){
		    totalRec++;
		    if(lineSp1[i] == '-'){

		    }else{
		     
			if(lineSp2[i] == '-'){
			 
			}else{
			    writtenRec++;

			    AlleleRecords arCurrent (false);
			    arCurrent.chri          = chrIdx;
			    arCurrent.coordinate    = coordC;
			    arCurrent.sizePops      = 1;
		
			    arCurrent.vectorAlleles = new vector<SingleAllele>  ();

			    //cout<<coordC<<"\t"<<lineSp1[i]<<lineSp2[i]<<endl;
			    //string toprintS=chrname+"\t"+stringify(coordC)+"\t";
			    SingleAllele sample;
			    char cRef= char(toupper(lineSp1[i]));
			    char cAlt= char(toupper(lineSp2[i]));
			    //string toprint;
			    if(cRef == cAlt){
				//toprintS+=stringify(cRef)+",N\t";
				//toprint = "1,0";
				sample.setRefCount(1); sample.setAltCount(0);

				arCurrent.ref           = cRef;
				arCurrent.alt           = 'N';				     
			    }else{
				//toprintS+=stringify(cRef)+","+stringify(cAlt)+"\t";
				//toprint = "0,1";
				sample.setRefCount(0); sample.setAltCount(1);
				arCurrent.ref           = cRef;
				arCurrent.alt           = cAlt;
			    }
			 

			 
			    bool cpgFlag=false;

			    if( (lastC+1) == coordC){
				if(lastCharWasC && (cRef == 'G' || cAlt == 'G' )){
				    cpgFlag=true;				 
				}
			    }
			    //sample.setIsCpg(cpgEPO);
			    arCurrent.vectorAlleles->push_back(root);
			    arCurrent.vectorAlleles->push_back(anc);
			    //arToWrite.vectorAlleles->push_back(sample);

			    if(cpgFlag){
				if(lastToPrintSet){
				    lastToPrint.vectorAlleles->at(2).setIsCpg(true);//setting CpG for previous record
				    if(!gw->writeAlleleRecord(&lastToPrint)){
					cerr<<"AXT2ACF: error writing header "<<endl;
					exit(1);
				    }
				}
				sample.setIsCpg(true);
				arCurrent.vectorAlleles->push_back(sample);
			     
				if(!gw->writeAlleleRecord(&arCurrent)){//printing for current record
				    cerr<<"AXT2ACF: error writing header "<<endl;
				    exit(1);
				}
				lastToPrintSet=false;//was flushed
				// if(!lastToPrintS.empty())
				// 	 cout<<lastToPrintS<<":1"<<endl; //print last line			     
				// toprintS+="0,0:0\t0,0:0\t"+toprint+":1"; //print current line
				// cout<<toprintS<<endl;                    //print current line
				// lastToPrintS="";
			    }else{

				if(lastToPrintSet){
				    lastToPrint.vectorAlleles->at(2).setIsCpg(false);//setting CpG for previous record
				    if(!gw->writeAlleleRecord(&lastToPrint)){
					cerr<<"AXT2ACF: error writing header "<<endl;
					exit(1);
				    }
				}

				sample.setIsCpg(false);
				arCurrent.vectorAlleles->push_back(sample);
				lastToPrint    = arCurrent;
				lastToPrintSet = true;

				// if(!lastToPrintS.empty())
				// 	 cout<<lastToPrintS<<":0"<<endl;   //print last line
				// toprintS+="0,0:0\t0,0:0\t"+toprint;   //store current line
				// lastToPrintS=toprintS;		       			 
			    }

			    //cout<<toPrintS<<endl;
			 
			    //cpg
			    if(cRef == 'C' || cAlt == 'C')
				lastCharWasC = true;
			    else
				lastCharWasC = false;

			    lastC = coordC;
			}
			coordC++;
		    }

		}//end loop

		if(lastToPrintSet){
		    lastToPrint.vectorAlleles->at(2).setIsCpg(false);//setting CpG for previous record
		 
		    if(!gw->writeAlleleRecord(&lastToPrint)){
			cerr<<"AXT2ACF: error writing header "<<endl;
			exit(1);
		    }
		}
		// if(!lastToPrintS.empty())
		// 	 cout<<lastToPrintS<<":0"<<endl;

		if((coordC-1) != endC){
		    cerr<<"The block did not end with the proper coordinate got:"<<coordC<<" expected "<<endC<<endl;
		    return 1;
		}

		// cout<<startC<<"\t"<<endC<<endl;
		// return 1;
		getline (myAXTFile,line);	 
	    }else{
		cerr<<"Wrong number of fields for line "<<line<<endl;
		return 1;
	    }
	}
	myAXTFile.close();
    }else{
	cerr << "Unable to open file "<<axtfile<<endl;
	return 1;
    }


    // string line;
    // bool hasData= getline (myfile,line);


    // while (hasData ){
    // 	totalRec++;

    // 	// while(btr.hasData()){
    // 	// BAMTableObj * toprint=btr.getData();
    // 	vector<string> token= allTokens(line,'\t');
    // 	// cout<<"line "<<line<<endl;
    // 	string chr    = token[1];
    // 	unsigned int pos = destringify<unsigned int>(token[2]);
    // 	string genotype    = token[3];

    // 	// if(passedFilters(toprint,filtersVCF)){

    // 	if(firstLine){
    // 	    firstLine=false;
    // 	    if(epoFileB){

    // 		rtEPO = new ReadTabix( epoFile.c_str()  , 
    // 				       epoFileidx.c_str()  , 
    // 				       chr,
    // 				       int(pos),INT_MAX ); //the destructor should be called automatically

    // 		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
    // 	    }
    // 	}

    // 		if(epoFileB){

    // 	if(!lineLeftEPO){
    // 	    cerr<<"Error, no data in the EPO file "<< chr <<":"<< int(pos) <<endl;
    // 	    return 1;
    // 	}

    // 	if(epoChr != chr){
    // 	    cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<chr<<endl;
    // 	    return 1;
    // 	}


    // 	// cout<<"2"<<endl;
    // 	while(epoCoord != pos){
    // 	    if(epoCoord > pos){
    // 		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(line)<<"\t"<<lineFromEPO<<endl;
    // 		return 1;
    // 	    }

    // 	    if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
    // 		rtEPO->repositionIterator(chr , int(pos),INT_MAX);
    // 	    }

	    
    // 	    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);

    // 	    // lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
    // 	    // vector<string> fieldsEPO = allTokens(lineFromEPO,'\t');
    // 	    // epoChr                   = fieldsEPO[0];
    // 	    // epoCoord                 = string2uint(fieldsEPO[1]);					
    // 	    // if(fieldsEPO[9] == "1")
    // 	    //     cpgEPO=true;		    
    // 	    // else
    // 	    //     cpgEPO=false;		    
    // 	    // if(ancAllele){
    // 	    //     allel_chimp = fieldsEPO[3][0];//inferred ancestor
    // 	    // }else{
    // 	    //     allel_chimp = fieldsEPO[4][0];//chimp;
    // 	    // }

    // 	    // if(!lineLeftEPO){
    // 	    //     cerr<<"Error, missing data in the EPO file"<<*toprint<<endl;
    // 	    //     return 1;
    // 	    // }
    // 	}


    // 	if(epoCoord != pos){
    // 	    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<line<<"\t"<<lineFromEPO<<endl;
    // 		return 1;
    // 	}
    // 		}
    // 	//int refIdx  =base2int(allel_ref);
    // 	// int chimpIdx=base2int(allel_chimp);
    // 	char alt='N';
    // 	string s="ACGT";
    // 	string chimpString;
    // 	string ancString;
    // 	SingleAllele root;
    // 	SingleAllele anc;

    // 	//need an allele that is A,C,G,T
    // 	//cout<<pos<<"\t"<<allel_ref<<endl;
    // 	if(!isResolvedDNA(allel_ref)){
    // 	    goto nextline;
    // 	}

    // 	// cout<<"3"<<endl;

    // 	//determine alternative allele
    // 	for(int i=0;i<4;i++){
    // 	    if(s[i] != allel_ref){
    // 		if(genotype.find(s[i]) != string::npos){ //new non-ref allele found
    // 		    //toprint->hasAllele(i+1) ){
    // 		    alt=s[i];
    // 		}
    // 	    }
    // 	}

    // 	// cout<<*toprint<<endl;
    // 	// cout<<lineFromEPO<<endl;
		
    // 	// SingleAllele root;
    // 	// SingleAllele anc;
    // 	if(!epoFileB){ //no epo file
    // 	    // chimpString="0,0:0";
    // 	    // ancString="0,0:0";
    // 	    root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
    // 	    anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

    // 	}else{
    // 	    //unresolved ancestral allele (A,C,G,T)
    // 	    if(!isResolvedDNA(allel_chimp)){
    // 		//chimpString="0,0:0"; 
    // 		root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
    // 	    }
    // 	    //resolved ancestral allele
    // 	    else{
    // 		if(allel_chimp == allel_ref){//no diff between chimp and reference
    // 		    //chimpString="1,0:"+string(cpgEPO?"1":"0");
    // 		    root.setRefCount(1); root.setAltCount(0);  root.setIsCpg(cpgEPO); 

    // 		}else{
    // 		    if(alt == 'N'){//no alt defined, the chimp becomes the alt			    
    // 			alt = allel_chimp;
    // 			//chimpString="0,1:"+string(cpgEPO?"1":"0");
    // 			root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 							
    // 		    }else{
    // 			if(alt == allel_chimp){//alt is chimp
    // 			    //chimpString="0,1:"+string(cpgEPO?"1":"0");
    // 			    root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 
    // 			}else{ //tri-allelic site, discard
    // 			    //continue;
    // 			    goto nextline;
    // 			}
    // 		    }
    // 		}
    // 	    }

    
    // 	    if(!isResolvedDNA(allel_anc)){
    // 		//ancString="0,0:0"; 
    // 		root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 
    // 	    }
    // 	    //resolved ancestral allele
    // 	    else{
    // 		if(allel_anc == allel_ref){//no diff between ancestor and reference
    // 		    //ancString="1,0:"+string(cpgEPO?"1":"0");
    // 		    anc.setRefCount(1); anc.setAltCount(0);  anc.setIsCpg(cpgEPO); 

    // 		}else{
    // 		    if(alt == 'N'){//no alt defined, the ancestor becomes the alt			    
    // 			alt = allel_anc;
    // 			//ancString="0,1:"+string(cpgEPO?"1":"0");			    
    // 			anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 

    // 		    }else{
    // 			if(alt == allel_anc){//alt is ancestor
    // 			    //ancString="0,1:"+string(cpgEPO?"1":"0");
    // 			    anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 
    // 			}else{ //tri-allelic site, discard
    // 			    //continue;
    // 			    goto nextline;
    // 			}
    // 		    }
    // 		}
    // 	    }
    // 	}

    // 	// cout<<"4"<<endl;





    // 	if( (alt!='N' && genotype.find(alt) != string::npos) )   // has alternative
    // 	    {//no alt in bam table

    // 	    int refCount=0;
    // 	    int altCount=0;
    // 	    string rr = stringify(allel_ref)+ stringify(allel_ref);
    // 	    string ra = stringify(allel_ref)+ stringify(alt);
    // 	    string ar = stringify(alt)      + stringify(allel_ref);
    // 	    string aa = stringify(alt)      + stringify(alt);

    // 	    if(genotype == rr ){
    // 	    	refCount=2;
    // 	    }else{
    // 	    	if(genotype == aa){
    // 	    	    altCount=2;
    // 	    	}else{
    // 	    	    if( (genotype == ra) ||
    // 	    		(genotype == ar) ){
    // 	    		refCount=1;
    // 	    		altCount=1;
    // 	    	    }else{
    // 	    		//goto nextline;
    // 			cerr<<"Error: Potential error in the 23 and me file where the reference allele is not there "<<line<<"\t"<<lineFromEPO<<endl; //error
    // 			goto nextline;
    // 	    	    }			
    // 	    	}
    // 	    }

    // 	    //cout<<refIdx<<endl;
    // 	    writtenRec++;


    // 	    AlleleRecords arToWrite (false);
    // 	    arToWrite.chri          = chr2index[chr];
    // 	    arToWrite.coordinate    = pos;
    // 	    arToWrite.sizePops      = 1;
    // 	    arToWrite.ref           = allel_ref;
    // 	    arToWrite.alt           = alt;		
		
    // 	    arToWrite.vectorAlleles = new vector<SingleAllele>  ();
    // 	    // SingleAllele root;
    // 	    // SingleAllele anc;
    // 	    SingleAllele sample (refCount, altCount, 0);//no CpG provided
		
    // 	    arToWrite.vectorAlleles->push_back(root);
    // 	    arToWrite.vectorAlleles->push_back(anc);
    // 	    arToWrite.vectorAlleles->push_back(sample);
	    
    // 	    if(!gw->writeAlleleRecord(&arToWrite)){
    // 		cerr<<"Vcf2ACF: error writing header "<<endl;
    // 		exit(1);
    // 	    }	    
    // 	    // cout<<chr<<"\t"<< pos<<"\t"<<
    // 	    // 	allel_ref<<","<<
    // 	    // 	alt<<"\t"<<
    // 	    // 	chimpString<<"\t"<<
    // 	    // 	ancString<<"\t"<<
    // 	    // 	refCount<<","<<
    // 	    // 	altCount<<
    // 	    // 	":"<<("0")<<endl;//no cpg info

    // 	}else{
    // 	    goto nextline;
    // 	    //continue;//triallelic, skip
    // 	}
	
    //     nextline:
        
    //     hasData = getline (myfile,line);
	
    // }

    delete(gw);
    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRec<<" records, wrote "<<writtenRec<<" terminated gracefully"<<endl;

    return 0;
}

