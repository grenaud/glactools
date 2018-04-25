
#include "BPLINK2ACF.h"


// #include "BAMTableObj.h"
// #include "BAMTABLEreader.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;
uint8_t mask4=192; //11000000
uint8_t mask3= 48; //00110000
uint8_t mask2= 12; //00001100
uint8_t mask1=  3; //00000011

uint8_t masks  [4] = {mask1,mask2,mask3,mask4};

void BPLINK2ACF::setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO){

    //lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
    lineLeftEPO=(rtEPO->readLineKS(  ));
    //cerr<<kstringPtrEPO->s<<endl;
    if(!lineLeftEPO){
	cerr<<"Error, missing data in the EPO file"<<endl;
	exit(1);
    }

    char *p, *q;
    int i;

    for (p = kstrtok(kstringPtrEPO->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
	q = (char*)aux.p;
	*q = 0;
	//cout<<i<<" >"<<p<<"<  "<<endl;
	if(i==0){//chr
	    epoChr                   = string(p);
	    continue;
	}
	if(i==1){//coord
	    epoCoord                 = (unsigned int)strtoul(p, NULL, 0);//strtoul(p);
	    continue;
	}
	if(i==2){//human ref
	    allel_ref   = p[0];
	    continue;
	}
	if(i==3){//ancestor
	    allel_anc     = p[0];
	    continue;
	}
	if(i==4){//chimp
	    allel_chimp   = p[0];
	    continue;
	}
	if(i<9) continue;
	if(i==9){
	    if(strcmp(p,"1")==0){
		cpgEPO=true;		    
	    }else{
		cpgEPO=false;
	    }
	    continue;
	}

	break;
    }

    // vector<string> fieldsEPO  = allTokens(lineFromEPO,'\t');
    // epoChr                   = fieldsEPO[0];
    // epoCoord                 = string2uint(fieldsEPO[1]);					
    // if(fieldsEPO[9] == "1")
    // 	cpgEPO=true;		    
    // else
    // 	cpgEPO=false;		    


    // allel_ref   = fieldsEPO[2][0];//reference allele
    // allel_anc   = fieldsEPO[3][0];//inferred ancestor
    // allel_chimp = fieldsEPO[4][0];//chimp;


}



BPLINK2ACF::BPLINK2ACF(){


}

BPLINK2ACF::~BPLINK2ACF(){

}


string BPLINK2ACF::usage() const{
    const string usage=string("")+ " bplink2acf [options] <BPLINK prefix> "+"\n"+

	                      "\nThis program convert BPLINK files into ACF (prints to the stdout)\n"+
			      "\nPlease note that since only variants are normally present in PLINK files, the CpG flag will not be set\n"+
"\nThe <BPLINK prefix> must include the following files:\n"+
"\t<BPLINK prefix>.bed\n"+
"\t<BPLINK prefix>.bim\n"+
"\t<BPLINK prefix>.fam\n"+
	"\n"+

	
			      "\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
			      "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+

			      "\t"+"--epo [EPO file]"       +"\t" +"Use file as EPO alignment to set the (default: none)\n"+   
			      "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   
		  
			      "";
    return usage;
}

int BPLINK2ACF::run(int argc, char *argv[]){

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


    if(lastOpt != (argc-1)){
	cerr<<"The last argument is <BPLINK prefix> "<<endl;
	return 1;		
    }


    if(fastaIndex.size()==0 ){
	cerr<<"Must specify fai file "<<endl;
	return 1;	
    }

    if( !strEndsWith(fastaIndex,".fai") ){
	cerr<<"fasta index must end with .fai file "<<endl;
	return 1;	
    }
    fastaFile = fastaIndex.substr(0,fastaIndex.size()-4);
    
    if(!epoFileB){
	if(!isFile( fastaFile ) ){
	    cerr<<"If you do not specify EPO file, we need the fasta index  as this file format does not have the reference information"<<endl<<"fasta file "<<fastaFile<<" does not exist"<<endl;
	    return 1;
	}
    	//cerr<<"Must specify EPO file as this file format does not have the reference information "<<endl;	
    	//return 1;	
    }

    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    readFastaIndex(fastaIndex,chrFound,genomeLength);

    //cerr<<"Filter used: "<<*filtersVCF<<endl;
    //VCFreader vcfr   (string(argv[argc-3]),5);
    // BAMTABLEreader btr  (),5);
    string inputFileEIGprefix = string(argv[lastOpt]);
    string bedFile   = string(argv[argc-1])+".bed";
    string bimFile   = string(argv[argc-1])+".bim";
    string famFile   = string(argv[argc-1])+".fam";
    // string genoFile   = inputFileEIGprefix+".geno";
    // string snpFile    = inputFileEIGprefix+".snp";
    // string indFile    = inputFileEIGprefix+".ind";

   if(!isFile(bedFile)){
        cerr<<".bed file "<<bedFile<<" does not exist"<<endl;
        exit(1);
    }

   if(!isFile(bimFile)){
        cerr<<".bim file "<<bimFile<<" does not exist"<<endl;
        exit(1);
   }
   
   if(!isFile(famFile)){
        cerr<<".fam file "<<famFile<<" does not exist"<<endl;
        exit(1);
    }

   igzstream myfileBED;
   myfileBED.open(bedFile.c_str(), ios::in | ios::binary);
   
   if (!myfileBED.good()){
       cerr << "Unable to open file "<<bedFile<<endl;
       return 1;
   }
   
   igzstream myfileBIM;
   myfileBIM.open(bimFile.c_str(), ios::in);
   
   if (!myfileBIM.good()){
       cerr << "Unable to open file "<<bimFile<<endl;
       return 1;
   }
   
   igzstream myfileFAM;
   myfileFAM.open(famFile.c_str(), ios::in);
   
   if (!myfileFAM.good()){
       cerr << "Unable to open file "<<famFile<<endl;
       return 1;
   }
    


    
    BamTools::Fasta fastaReference;
   
    string epoFileidx = epoFile+".tbi";

    string epoChr;
    unsigned int epoCoord;

    ReadTabix * rtEPO =NULL;
    //string lineFromEPO;
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
    header+="#BPLINK2ACF:\n";
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



    header+="#chr\tcoord\tREF,ALT\troot\tanc";
    bool foundPop=false;
    bool foundAnc=false;
    int  indexAnc=-1;
    bool foundRoot=false;
    int  indexRoot=-1;

    string lineFAM;
    int sizePopRead=0;
    //int allPopFound=0;
    
    int indexPop=-1;
    while ( getline (myfileFAM,lineFAM)){   
	//header<<line<<endl;
	trimWhiteSpacesBothEnds(&lineFAM);
	vector<string> fields = allTokensWhiteSpaces( lineFAM );
	if(fields.size() != 6){
	    cerr<<"Line "<<lineFAM<<" should have 6 TAB/space deliminated fields, found "<<fields.size()<<" tab/space deliminated fields"<<endl;
	    return 1;
	}
	//allPopFound++;
	indexPop++;
	if(fields[0] == "root"){ foundRoot=true; indexRoot=indexPop; continue; }
	if(fields[0] == "anc") { foundAnc=true;  indexAnc =indexPop; continue; }
       
	sizePopRead++;
	header+="\t"+fields[1];
	foundPop=true;
    }
    myfileFAM.close();

    if( (foundAnc || foundRoot) && epoFileB){
	cerr<<"ERROR: we found a pop called root and anc in the FAM file:"<<famFile<<" but the EPO file is used. We cannot determine from where to take the ancestral/root information from, exiting"<<endl;
	return 1;


    }
    if(!foundPop){
	cerr<<"ERROR: did not find any data in the individual file"<<endl;
	return 1;
    }

    header+="\n";

    int bytesToRead=ceil(  double(sizePopRead)/4.0 );


    //cerr<<header<<endl;
    //return 1;
    
    GlacWriter * gw=NULL;
    //cerr<<"sizePopRead "<<sizePopRead<<endl;
    
    gw = new GlacWriter(sizePopRead,    //gp.getSizePops(),
                        false,          //gp.isGLFormat(),
                        2,              //gp.isACFormat()?2:1,
			1,              //compression threads
                        uncompressed);
    
    
    if(!gw->writeHeader(header)){
	cerr<<"GlacViewer: error writing header "<<endl;
	exit(1);
    }

    //checking magic numbers
    char c;
    myfileBED.read((char*)&c,sizeof(char));
    if(int(c) != 108){
	cerr<<"BPLINK2ACF: error in magic numbers found "<<int(c)<<" expected 108"<<endl;
	exit(1);	
    }
    myfileBED.read((char*)&c,sizeof(char));
    if(int(c) != 27){
	cerr<<"BPLINK2ACF: error in magic numbers found "<<int(c)<<" expected 27"<<endl;
	exit(1);	
    }
    myfileBED.read((char*)&c,sizeof(char));
    if(int(c) != 1){
	cerr<<"BPLINK2ACF: error in magic numbers found "<<int(c)<<" expected 1"<<endl;
	exit(1);	
    }



    
    string lineBIM;
    //
    while ( getline (myfileBIM,lineBIM)){   
	//header<<line<<endl;
	trimWhiteSpacesBothEnds(&lineBIM);
	totalRec++;

	//write code to retrieve data from bed file here
    	vector<string> fields = allTokensWhiteSpaces( lineBIM );
	if(fields.size() != 6){
	    cerr<<"Line "<<lineBIM<<" should have 6 TAB/space deliminated fields, found "<<fields.size()<<" tab/space deliminated fields"<<endl;
	    return 1;
	}
	//token[0];//SNPname
	
    	string chr          = fields[0];//chr
	//1 - ID
	//2 pos in CM
    	//token[2];//genetic position
    	unsigned int pos    = destringify<unsigned int>(fields[3]);//physical position
    	char genotypeMin    = fields[4][0];//minor first
    	char genotypeMaj    = fields[5][0];//major second

	

    	if(firstLine){
    	    firstLine=false;
    	    if(epoFileB){

    		rtEPO = new ReadTabix( epoFile.c_str()  , 
    				       epoFileidx.c_str()  , 
    				       chr,
    				       int(pos),INT_MAX ); //the destructor should be called automatically
    		kstringPtrEPO = rtEPO->getKstringPtr();
    		memset(&aux, 0, sizeof(ks_tokaux_t));
    		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
    	    }else{
    		if ( !fastaReference.Open(fastaFile , fastaIndex) ){ 
    		    cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " <<fastaIndex<<endl;
    		    return 1;
    		}		
    	    }
    	}

    	if(epoFileB){
    	    if(!lineLeftEPO){
    		cerr<<"Error, no data in the EPO file "<< chr <<":"<< int(pos) <<endl;
    		return 1;
    	    }
	    
    	    if(epoChr != chr){
    		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
    		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
    		}	       	    
    		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
    	    }

    	    while(epoCoord != pos){
    		if(epoCoord > pos){
    		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(lineBIM)<<"\tEPO="<<kstringPtrEPO->s<<endl;
    		    return 1;
    		}		
    		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
    		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
    		}			    
    		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
    	    }

    	    if(epoCoord != pos){
    		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<lineBIM<<"\tEPO="<<kstringPtrEPO->s<<endl;
    		return 1;
    	    }
    	}else{ //end if epoFileB
    	    //allel_ref =
    	    if ( !fastaReference.GetBase(chr2index[chr], pos-1, allel_ref ) ) {
    		cerr << "glactools ERROR:  could not read reference base from FASTA file at chr:"<<chr2index[chr]<<" position:"<<(pos) << endl;
    		 return 1;
    	     }

    	}
	
    	char alt='N';
    	string s="ACGT";
    	string chimpString;
    	string ancString;
    	SingleAllele root;
    	SingleAllele anc;
    	AlleleRecords arToWrite (false);
    	//bool refIsAltfromSNPfile=false;	
  	char genotypeRef;
    	char genotypeAlt;
	bool majorIsRef=false;
	bool minorIsRef=false;
	int indexOfIndividual=0;	
	vector<SingleAllele> allexeptRootAnc;

	//need an allele that is A,C,G,T
    	if(!isResolvedDNA(genotypeMin)){
    	    goto nextline;
    	}
    	if(!isResolvedDNA(genotypeMaj)){
    	    goto nextline;
    	}
	//cerr<<"pos:"<<pos<<" ref:"<<allel_ref<<" maj:"<<genotypeMaj<<" min:"<<genotypeMin<<endl;
	if(isResolvedDNA(allel_ref)){

	    //need to set genotypeAlt
	    if(allel_ref == genotypeMaj){//if the reference is the Max
		majorIsRef=true;
		genotypeRef=allel_ref;
		if(isResolvedDNA(genotypeMin)){//minor is resolved
		    if(genotypeMaj!=genotypeMin){
			genotypeAlt = genotypeMin;
		    }else{
			genotypeAlt = 'N';
		    }			
		}else{//minor is not resolved, should not happen
		    genotypeAlt = 'N';
		}
		
		
	    }else{
		
		if(allel_ref == genotypeMin){//if the reference is the Max
		    minorIsRef=true;
		    genotypeRef=allel_ref;
		    
		    if(isResolvedDNA(genotypeMaj)){//major is resolved
			if(genotypeMaj!=genotypeMin){//difference
			    genotypeAlt = genotypeMaj;
			}else{
			    genotypeAlt = 'N';//should not happen
			}			
		    }else{//minor is not resolved, should not happen
			genotypeAlt = 'N';
		    }
		    
		}else{
		    //ref not found
		    cerr<<"bplink2acf: WARNING: The reference allele between the EPO/FASTA ("<<allel_ref<<") was not found in the the .bim file (found: "<<genotypeMaj<<","<<genotypeMin<<") at chr:pos "<<chr<<":"<<pos<<endl;
		    return 1;
		    goto nextline;
		}
		
	    }
	    
	}else{//if not resolved for allele refs
	    goto nextline;
	}
	//cerr<<"pos:"<<pos<<" majorIsRef:"<<majorIsRef<<" minorIsRef:"<<minorIsRef<<endl;
	
	if(!majorIsRef && !minorIsRef){
	    cerr<<"bplink2acf: WARNING: The reference allele between the EPO/FASTA ("<<allel_ref<<") was not found in the the .bim file (found: "<<genotypeMaj<<","<<genotypeMin<<") at chr:pos "<<chr<<":"<<pos<<endl;
	    goto nextline;
	}



	if(!epoFileB){ //no epo file

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

	
    	if(alt!='N'){
    	    if(genotypeAlt != alt ){   // has alternative
    		cerr<<"bplink2acf: WARNING The alternative allele between the EPO and the .bim file disagrees, could be triallelic at chr:pos "<<chr<<":"<<thousandSeparator(pos)<<", skipping"<<endl;	
    		goto nextline;
    	    }	    
    	}



	
    	arToWrite.chri          = chr2index[chr];
    	arToWrite.coordinate    = pos;
    	arToWrite.sizePops      = uint32_t(sizePopRead);
    	arToWrite.ref           = genotypeRef;
    	arToWrite.alt           = genotypeAlt;
		
    	arToWrite.vectorAlleles = new vector<SingleAllele>  ();
    	// SingleAllele root;
    	// SingleAllele anc;
	//cerr<<"bytesToRead "<<bytesToRead<<endl;

    	for(int i=0;i<bytesToRead;i++){
	    unsigned char data;
	    myfileBED.read((char*)&data,sizeof(char));
	    //mask and shift
	    //cerr<<"data: "<<var2binary(data)<<" "<<indexOfIndividual<<" "<<(indexOfIndividual%4)<<" "<<sizePopRead<<endl;
	    for(int indBytes=0;indBytes<4;indBytes++){
		if(indexOfIndividual>indexPop)
		    break;
		unsigned char dataInd;
		//cerr<<"mask: "<<var2binary(masks[ (indexOfIndividual%4) ])<<endl;
		dataInd = (data & masks[ (indexOfIndividual%4) ]);
		// cerr<<"dataInd1: "<<var2binary(dataInd)<<endl;
		// cerr<<indexOfIndividual<<" "<<(indexOfIndividual%4)<<" "<<((indexOfIndividual%4)*2)<<endl;
		dataInd = (dataInd >> ( (indexOfIndividual%4)*2  ));
		//cerr<<"dataInd2: "<<var2binary(dataInd)<<endl;
	    
		int g =int(dataInd);
		//cerr<<"g: "<<g<<endl;

		int refCount=0;
		int altCount=0;
	    
		//homo for minor
		if(g == 0){
		    if(majorIsRef){
			refCount = 0;  altCount = 2;
		    }else{//minor is ref
			refCount = 2;  altCount = 0;
		    }
		}
		
		//missing
		if(g == 1){
		    refCount = 0;  altCount = 0;
		}
	    
		//het
		if(g == 2){
		    refCount = 1;  altCount = 1;//refIsAltfromSNPfile does not matter
		}
	    
		//homo minor
		if(g == 3){		
		    if(majorIsRef){
			refCount = 2;  altCount = 0;
		    }else{//minor is ref
			refCount = 0;  altCount = 2;
		    }
		}

		if( indexOfIndividual ==  indexAnc){
		    anc.setRefCount(refCount);   anc.setAltCount(altCount);  anc.setIsCpg(false); 
		}else{
    
		    if( indexOfIndividual ==  indexRoot){
			root.setRefCount(refCount);  root.setAltCount(altCount); root.setIsCpg(false); 
		    }else{

			SingleAllele sample (refCount, altCount, 0);//no CpG provided
			allexeptRootAnc.push_back(sample);

		    }
		}
		//cerr<<"indexOfIndividual "<<indexAnc<<" "<<indexRoot<<" "<<indexOfIndividual<<" "<< indexPop <<" sizes:"<<allexeptRootAnc.size()<<endl;
		indexOfIndividual++;
	    }
    	}//for each char in lineG


	//cerr<<"writing sizes:"<<allexeptRootAnc.size()<<endl;
    	arToWrite.vectorAlleles->push_back(root);
    	arToWrite.vectorAlleles->push_back(anc);
	for(unsigned int j=0;j<allexeptRootAnc.size();j++){
	    //cerr<<j<<" "<<allexeptRootAnc[j]<<endl;
	    arToWrite.vectorAlleles->push_back(allexeptRootAnc[j]);
	}
	
    	if(!gw->writeAlleleRecord(&arToWrite)){
    	    cerr<<"bplink2acf: error writing header "<<endl;
    	    exit(1);
    	}	    
	
    	writtenRec++;
	
    nextline:
	continue;
    }
    myfileBIM.close();
    myfileBED.close();



    
    if(rtEPO != NULL)
	delete(rtEPO);
    
    delete(gw);

    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRec<<" records, wrote "<<writtenRec<<" terminated gracefully"<<endl;

    return 0;
}











    /*





    
    string lineG;//geno
    string lineS;//snp
    
    bool hasDataG= (bool)getline (myfileG,lineG);
    bool hasDataS= (bool)getline (myfileS,lineS);

    while (hasDataG ){
    	if(!hasDataS){	    
    	    cerr<<"bplink2acf: data seems to be missing in the .snp file, stopped at line in .geno "<<lineG<<endl;
    	    return 1;
    	}


	
    	totalRec++;
    	// cout<<lineS<<endl;
    	// cout<<lineG<<endl;
    	trimWhiteSpacesBothEnds(&lineS);
	
    	vector<string> token= allTokensWhiteSpaces(lineS);
	
    	if(token.size() != 6){
    	    cerr<<"Error, SNP line should contain 6 fields"<<endl;
    	    return 1;
    	}
	
    	//token[0];//SNPname
    	string chr          = token[1];//chr
    	//token[2];//genetic position
    	unsigned int pos    = destringify<unsigned int>(token[3]);//physical position
    	char genotypeRef    = token[4][0];
    	char genotypeAlt    = token[5][0];

	

    	if(firstLine){
    	    firstLine=false;
    	    if(epoFileB){

    		rtEPO = new ReadTabix( epoFile.c_str()  , 
    				       epoFileidx.c_str()  , 
    				       chr,
    				       int(pos),INT_MAX ); //the destructor should be called automatically
    		kstringPtrEPO = rtEPO->getKstringPtr();
    		memset(&aux, 0, sizeof(ks_tokaux_t));
    		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
    	    }else{
    		if ( !fastaReference.Open(fastaFile , fastaIndex) ){ 
    		    cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " <<fastaIndex<<endl;
    		    return 1;
    		}		
    	    }
    	}

    	if(epoFileB){

    	    if(!lineLeftEPO){
    		cerr<<"Error, no data in the EPO file "<< chr <<":"<< int(pos) <<endl;
    		return 1;
    	    }

    	    if(epoChr != chr){
    		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
    		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
    		}	       	    
    		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
    	    }



    	    while(epoCoord != pos){
    		if(epoCoord > pos){
    		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(lineS)<<"\tEPO="<<kstringPtrEPO->s<<endl;
    		    return 1;
    		}
		
    		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
    		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
    		}
		
	    
    		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
    	    }


    	    if(epoCoord != pos){
    		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<lineS<<"\tEPO="<<kstringPtrEPO->s<<endl;
    		return 1;
    	    }
    	}else{ //end if epoFileB

    	    //allel_ref =
    	    if ( !fastaReference.GetBase(chr2index[chr], pos, allel_ref ) ) {
    		cerr << "glactools ERROR:  could not read reference base from FASTA file at chr:"<<chr2index[chr]<<" position:"<<(pos) << endl;
    		 return 1;
    	     }

    	}
	
    	char alt='N';
    	string s="ACGT";
    	string chimpString;
    	string ancString;
    	SingleAllele root;
    	SingleAllele anc;
    	AlleleRecords arToWrite (false);
    	bool refIsAltfromSNPfile=false;	
    	//need an allele that is A,C,G,T
    	if(!isResolvedDNA(genotypeRef)){
    	    goto nextline;
    	}
    	if(!isResolvedDNA(genotypeAlt)){
    	    goto nextline;
    	}

    	if(epoFileB){
    	    if(isResolvedDNA(allel_ref)){
    		if(allel_ref != genotypeRef){
    		    //cerr<<"eigen2acf: WARNING: The reference allele between the EPO ("<<allel_ref<<") and the ref from the .geno file ("<<genotypeRef<<") disagrees at chr:pos "<<chr<<":"<<pos<<endl;
    		    if(allel_ref != genotypeAlt){
    			cerr<<"bplink2acf: WARNING: The reference allele between the EPO ("<<allel_ref<<") and the alt from the .geno file ("<<genotypeRef<<") also disagrees at chr:pos "<<chr<<":"<<pos<<", skiping"<<endl;			
    			goto nextline;
    		    }else{
    			refIsAltfromSNPfile=true;
    			char c_     = genotypeRef;
    			genotypeRef = genotypeAlt;
    			genotypeAlt = c_;
    		    }
		    
    		}else{
    		    //fine
		    
    		}
    	    }else{
    		goto nextline;
    	    }
	    
    	}else{
	    if(isResolvedDNA(allel_ref)){
    		if(allel_ref != genotypeRef){
    		    //cerr<<"eigen2acf: WARNING: The reference allele between the EPO ("<<allel_ref<<") and the ref from the .geno file ("<<genotypeRef<<") disagrees at chr:pos "<<chr<<":"<<pos<<endl;
    		    if(allel_ref != genotypeAlt){
    			cerr<<"bplink2acf: WARNING: The reference allele between the EPO ("<<allel_ref<<") and the alt from the .geno file ("<<genotypeRef<<") also disagrees at chr:pos "<<chr<<":"<<pos<<", skiping"<<endl;			
    			goto nextline;
    		    }else{
    			refIsAltfromSNPfile=true;
    			char c_     = genotypeRef;
    			genotypeRef = genotypeAlt;
    			genotypeAlt = c_;
    		    }
		    
    		}else{
    		    //fine
		    
    		}
    	    }else{
    		goto nextline;
    	    }
    	}
	   
    	//todo

    // 	    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
    // 	cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file at chr "<<pileupData.RefId<<" position "<<(posAlign-1) << endl;
    // 	exit(1);
    // }

    	//determine alternative allele
    	// for(int i=0;i<4;i++){
    	//     if(s[i] != allel_ref){
    	// 	if(genotype.find(s[i]) != string::npos){ //new non-ref allele found
    	// 	    //toprint->hasAllele(i+1) ){
    	// 	    alt=s[i];
    	// 	}
    	//     }
    	// }


    	if(!epoFileB){ //no epo file

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

	
    	if(alt!='N'){
    	    if(genotypeAlt != alt ){   // has alternative
    		cerr<<"bplink2acf: WARNING The alternative allele between the EPO and the .geno file disagrees, could be triallelic at chr:pos "<<chr<<":"<<thousandSeparator(pos)<<", skipping"<<endl;		
    		goto nextline;
    	    }	    
    	}





    	arToWrite.chri          = chr2index[chr];
    	arToWrite.coordinate    = pos;
    	arToWrite.sizePops      = uint32_t(sizePopRead);
    	arToWrite.ref           = genotypeRef;
    	arToWrite.alt           = genotypeAlt;
		
    	arToWrite.vectorAlleles = new vector<SingleAllele>  ();
    	// SingleAllele root;
    	// SingleAllele anc;
	
    	arToWrite.vectorAlleles->push_back(root);
    	arToWrite.vectorAlleles->push_back(anc);

    	for(unsigned int i=0;i<lineG.size();i++){
    	    char g = lineG[i];
    	    int refCount=0;
    	    int altCount=0;

	    
    	    //homo ref
    	    if(g == '2'){		
    		if(!refIsAltfromSNPfile){
    		    refCount = 2;  altCount = 0;
    		}else{
    		    refCount = 0;  altCount = 2;
    		}
    	    }

    	    //het
    	    if(g == '1'){
    		refCount = 1;  altCount = 1;//refIsAltfromSNPfile does not matter
    	    }
	    
    	    //homo alt
    	    if(g == '0'){
    		if(!refIsAltfromSNPfile){
    		    refCount = 0;  altCount = 2;
    		}else{
    		    refCount = 2;  altCount = 0;
    		}
    	    }

    	    if(g == '9'){
    		refCount = 0;  altCount = 0;
    	    }

	    
    	    SingleAllele sample (refCount, altCount, 0);//no CpG provided

    	    arToWrite.vectorAlleles->push_back(sample);
    	}//for each char in lineG
	
    	if(!gw->writeAlleleRecord(&arToWrite)){
    	    cerr<<"bplink2acf: error writing header "<<endl;
    	    exit(1);
    	}	    
	
    	writtenRec++;
	
    nextline:
        
        hasDataG = (bool)getline (myfileG,lineG);
        hasDataS = (bool)getline (myfileS,lineS);
	
    }//while hasData
*/
