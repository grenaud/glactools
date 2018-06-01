
#include "BEAGLE2GLF.h"


static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;

//returns PHRED of genotype likelihood but capped at 2000 (does not matter as they will be capped at 255)
inline int BEAGLE2GLF::phredcapped(const double p) const{
    if(p==0) return 2000;
    int t=int(-10.0*log10( p ));
    if(t>2000) return 2000;
    return t;
}


void BEAGLE2GLF::setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO){

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



BEAGLE2GLF::BEAGLE2GLF(){


}

BEAGLE2GLF::~BEAGLE2GLF(){

}


string BEAGLE2GLF::usage() const{
    const string usage=string("")+ " beagle2glf [options] <BEAGLE prefix> "+"\n"+
			      "\nThis program convert BEAGLE files into GLF (prints to the stdout)\n"+

	"\nThe <BEAGLE prefix> must include the following files:\n"+
	//"\t<BEAGLE prefix>.geno\n"+
"\t<BEAGLE prefix>.beagle.gz\n"+
"\t<BEAGLE prefix>.pos.gz\n"+
"\n"+
	
			      "\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
			      "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+

			      "\t"+"--epo [EPO file]"       +"\t" +"Use file as EPO alignment to set the (default: none)\n"+   
			      "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   
		  
			      "";
    return usage;
}

int BEAGLE2GLF::run(int argc, char *argv[]){

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
	cerr<<"The last argument is <BEAGLE prefix> "<<endl;
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

    //string genoFile   = inputFileEIGprefix+".geno";
    string beagleFile    = inputFileEIGprefix+".beagle.gz";
    string posFile       = inputFileEIGprefix+".pos.gz";

    
    string lineG;//geno    
    string lineB;//beagle
    string lineP;//pos

   // if(!isFile(genoFile)){
   //      cerr<<".geno file "<<genoFile<<" does not exist"<<endl;
   //      exit(1);
   //  }

   if(!isFile(beagleFile)){
        cerr<<".beagle.gz file "<<beagleFile<<" does not exist"<<endl;
        exit(1);
   }
   
   if(!isFile(posFile)){
        cerr<<".pos.gz file "<<posFile<<" does not exist"<<endl;
        exit(1);
    }

   // igzstream myfileG;
   // myfileG.open(genoFile.c_str(), ios::in);

   //  if (!myfileG.good()){
   //      cerr << "Unable to open file "<<genoFile<<endl;
   //      return 1;
   //  }

    igzstream myfileB;
    myfileB.open(beagleFile.c_str(), ios::in);

    if (!myfileB.good()){
        cerr << "Unable to open file "<<beagleFile<<endl;
        return 1;
    }

    igzstream myfileP;
    myfileP.open(posFile.c_str(), ios::in);

    if (!myfileP.good()){
        cerr << "Unable to open file "<<posFile<<endl;
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

    header+="#GLF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    header+="#PG:"+programLine+"\n";;
    header+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";;

    header+="#DATE: "+getDateString()+"\n";;
    header+="#BEAGLE2GLF:\n";
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


    int sizePopRead=0;
    header+="#chr\tcoord\tREF,ALT\troot\tanc";
    bool foundPop=false;
    //get header
    if(getline (myfileB,lineB)){
	vector<string> fields = allTokensWhiteSpaces( lineB );
	if(fields.size() < 3){
	    cerr<<"Line "<<lineB<<" should have 3 TAB/space deliminated fields, found "<<fields.size()<<" tab/space deliminated fields"<<endl;
	    return 1;
	}
	
	for(unsigned int p=3;p<fields.size();p+=3){
	    header+="\t"+fields[p];
	    if(fields[p] != fields[p+1] ){
		cerr<<"Line "<<lineB<<" should have 3 identical IDs for the pops"<<endl;
		return 1;
	    }
	    if(fields[p] != fields[p+2] ){
		cerr<<"Line "<<lineB<<" should have 3 identical IDs for the pops"<<endl;
		return 1;
	    }
	    sizePopRead++;
	}

	foundPop=true;
    }else{
	cerr<<"Did not find any data in the BEAGLE file"<<endl;
	return 1;
    }
    // cerr<<sizePopRead<<endl;
    // return 1;

    //get header
    if(getline (myfileP,lineP)){
    }else{
	cerr<<"Did not find any data in the BEAGLE file"<<endl;
	return 1;
    }

    // string linei;
    // int sizePopRead=0;
    // while ( getline (myfileI,linei)){   
    // 	//header<<line<<endl;
    // 	trimWhiteSpacesBothEnds(&linei);
    // 	vector<string> fields = allTokensWhiteSpaces( linei );
    // 	if(fields.size()!=3){
    // 	    cerr<<"Line "<<linei<<" should have 3 TAB/space deliminated fields, found "<<fields.size()<<" tab/space deliminated fields"<<endl;
    // 	    return 1;
    // 	}
    // 	sizePopRead++;
    // 	header+="\t"+fields[0];
    // 	foundPop=true;
    // }
    // myfileI.close();
    
    if(!foundPop){
	cerr<<"Did not find any data in the BEAGLE file"<<endl;
	return 1;
    }

    header+="\n";
    // cout<<header<<endl;
    // return 1;

    GlacWriter * gw=NULL;

    gw = new GlacWriter(sizePopRead,    //gp.getSizePops(),
                        true, //gp.isGLFormat(),
                        1,//gp.isACFormat()?2:1,
			1,//compression threads
                        uncompressed);
    

    if(!gw->writeHeader(header)){
	cerr<<"GlacViewer: error writing header "<<endl;
	exit(1);
    }



    


    
    bool hasDataB= (bool)getline (myfileB,lineB);
    bool hasDataP= (bool)getline (myfileP,lineP);

    while( hasDataP ){
	if(!hasDataB){	    
	    cerr<<"beagle2glf: data seems to be missing in the .beagle file, stopped at line in .pos "<<lineP<<endl;
	    return 1;
	}

	totalRec++;
	// cerr<<lineB<<endl;
	// cerr<<lineP<<endl;
	
	trimWhiteSpacesBothEnds(&lineB);
	trimWhiteSpacesBothEnds(&lineP);
	
	vector<string> tokenB = allTokensWhiteSpaces(lineB);
	vector<string> tokenP = allTokensWhiteSpaces(lineP);
	
	if(tokenP.size() < 2){
	    cerr<<"Error, POS line should contain 3 fields"<<lineP<<endl;
	    return 1;
	}
	
	//token[0];//SNPname
	string chr          = tokenP[0];//chr	
	//token[2];//genetic position
	unsigned int pos    = destringify<unsigned int>(tokenP[1]);//physical position
	int i_genotypeRef   = destringify<int>(tokenB[1]);
	int i_genotypeAlt   = destringify<int>(tokenB[2]);
	if( i_genotypeRef<0 ||  i_genotypeRef>3){
	    cerr<<"Error, POS line "<<lineP<<" contains an unknown genome"<<endl;
	    return 1;
	}
	if( i_genotypeAlt<0 ||  i_genotypeAlt>3){
	    cerr<<"Error, POS line "<<lineP<<" contains an unknown genome"<<endl;
	    return 1;
	}
	   

	char genotypeRef    = "ACGT"[i_genotypeRef];
	char genotypeAlt    = "ACGT"[i_genotypeAlt];

	
	
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
		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(lineP)<<"\tEPO="<<kstringPtrEPO->s<<endl;
		    return 1;
		}
		
		if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO->repositionIterator(chr , int(pos),INT_MAX);
		}
		
	    
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO);
	    }


	    if(epoCoord != pos){
		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<lineP<<"\tEPO="<<kstringPtrEPO->s<<endl;
		return 1;
	    }
	}else{ //end if epoFileB

	    //allel_ref =
	    if ( !fastaReference.GetBase(chr2index[chr], pos-1, allel_ref ) ) {
		cerr << "glactools ERROR:  could not read reference base from FASTA file at chr:"<<chr2index[chr]<<" position:"<<(pos) << endl;
		 return 1;
	     }

	}
	
	//to be safe
	allel_ref = toupper( allel_ref );
	allel_anc = toupper( allel_anc );
	


	
	char alt='N';
	string s="ACGT";
	string chimpString;
	string ancString;
	SingleGL root;
	SingleGL anc;
	AlleleRecords arToWrite (false);
	bool refIsAltfromSNPfile=false;	
	//need an allele that is A,C,G,T
	if(!isResolvedDNA(genotypeRef)){
	    goto nextline;
	}
	if(!isResolvedDNA(genotypeAlt)){
	    goto nextline;
	}

	if(isResolvedDNA(allel_ref)){
	    if(allel_ref != genotypeRef){
		//cerr<<"eigen2acf: WARNING: The reference allele between the EPO ("<<allel_ref<<") and the ref from the .geno file ("<<genotypeRef<<") disagrees at chr:pos "<<chr<<":"<<pos<<endl;
		if(allel_ref != genotypeAlt){
		    cerr<<"eigen2acf: WARNING: The reference allele between the EPO/FASTA ("<<allel_ref<<") and the alt from the .geno file ("<<genotypeRef<<") also disagrees at chr:pos "<<chr<<":"<<pos<<", skiping"<<endl;			
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
	
	//cerr<< genotypeRef<<" "<<genotypeAlt<<" "<<allel_ref<<" "<<allel_anc<<" "<<refIsAltfromSNPfile<<endl;


	if(!epoFileB){ //no epo file

	    // root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
	    // anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

	    root.setrrGL(0); root.setraGL(0);  root.setaaGL(0); root.setIsCpg(false);
	    anc.setrrGL(0);  anc.setraGL(0);   anc.setaaGL(0);  anc.setIsCpg(false);

	}else{
	    //unresolved ancestral allele (A,C,G,T)
	    if(!isResolvedDNA(allel_chimp)){
		//chimpString="0,0:0"; 
		//root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false);
		root.setrrGL(0); root.setraGL(0);  root.setaaGL(0); root.setIsCpg(false);
	    }
	    //resolved ancestral allele
	    else{
		if(allel_chimp == allel_ref){//no diff between chimp and reference
		    //chimpString="1,0:"+string(cpgEPO?"1":"0");
		    // root.setRefCount(1); root.setAltCount(0);  root.setIsCpg(cpgEPO); 
		    root.setrrGL(0); root.setraGL(255); root.setaaGL(255);  root.setIsCpg(cpgEPO);
		}else{
		    if(alt == 'N'){//no alt defined, the chimp becomes the alt			    
			alt = allel_chimp;
			//chimpString="0,1:"+string(cpgEPO?"1":"0");
			// root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO);
			root.setrrGL(255); root.setraGL(255); root.setaaGL(0);  root.setIsCpg(cpgEPO);
		    }else{
			if(alt == allel_chimp){//alt is chimp
			    //chimpString="0,1:"+string(cpgEPO?"1":"0");
			    //root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO);
			    root.setrrGL(255); root.setraGL(255); root.setaaGL(0);  root.setIsCpg(cpgEPO);
			}else{ //tri-allelic site, discard
			    //continue;
			    goto nextline;
			}
		    }
		}
	    }


	    if(!isResolvedDNA(allel_anc)){
		//ancString="0,0:0"; 
		//anc.setRefCount(0); anc.setAltCount(0);  anc.setIsCpg(cpgEPO); 
		anc.setrrGL(0);     anc.setraGL(0);      anc.setaaGL(0); anc.setIsCpg(false);

	    }
	    //resolved ancestral allele
	    else{
		if(allel_anc == allel_ref){//no diff between ancestor and reference
		    //ancString="1,0:"+string(cpgEPO?"1":"0");
		    //anc.setRefCount(1); anc.setAltCount(0);  anc.setIsCpg(cpgEPO); 
		    anc.setrrGL(0);     anc.setraGL(255);      anc.setaaGL(255); anc.setIsCpg(cpgEPO);

		}else{
		    if(alt == 'N'){//no alt defined, the ancestor becomes the alt			    
			alt = allel_anc;
			//ancString="0,1:"+string(cpgEPO?"1":"0");			    
			//anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 
			anc.setrrGL(255);     anc.setraGL(255);      anc.setaaGL(0); anc.setIsCpg(cpgEPO);

		    }else{
			if(alt == allel_anc){//alt is ancestor
			    //ancString="0,1:"+string(cpgEPO?"1":"0");
			    //anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO);
			    anc.setrrGL(255);     anc.setraGL(255);      anc.setaaGL(0); anc.setIsCpg(cpgEPO);

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
		cerr<<"beagle2glf: WARNING The alternative allele between the EPO and the .geno file disagrees, could be triallelic at chr:pos "<<chr<<":"<<thousandSeparator(pos)<<", skipping"<<endl;		
		goto nextline;
	    }	    
	}
	//cerr<< genotypeRef<<" "<<genotypeAlt<<" "<<allel_ref<<" "<<allel_anc<<" "<<refIsAltfromSNPfile<<endl;
	
	//	cerr<<genotypeRef<<"\t"<<genotypeAlt<<endl;


	


	arToWrite.chri          = chr2index[chr];
	arToWrite.coordinate    = pos;
	arToWrite.sizePops      = uint32_t(sizePopRead);
	arToWrite.ref           = genotypeRef;
	arToWrite.alt           = genotypeAlt;

	
	//todo
	arToWrite.vectorGLs = new vector<SingleGL>  ();
	//arToWrite.vectorAlleles = new vector<SingleAllele>  ();
	// SingleAllele root;
	// SingleAllele anc;
	
	arToWrite.vectorGLs->push_back(root);
	arToWrite.vectorGLs->push_back(anc);
	//https://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf
	
	for(unsigned int i=3;i<tokenB.size();i+=3){
	    // cerr<<i<<" "<<tokenB.size()<<" "<<arToWrite.vectorGLs->size()<<endl;
	    // cerr<<(tokenB[i+0])<<endl;
	    // cerr<<(tokenB[i+1])<<endl;
	    // cerr<<(tokenB[i+2])<<endl;

	    double g_AA = destringify<double>(tokenB[i+0]);
	    double g_AB = destringify<double>(tokenB[i+1]);
	    double g_BB = destringify<double>(tokenB[i+2]);
	    
	    double sum=g_AA+g_AB+g_BB;
	    if(sum<0.95 || sum>1.05){
		cerr<<"BEAGLE2GLF ERROR: the genotype likelihood for individual "<<(i/3)<<" do not sum up to ~1 at "<<lineB<<"\t sum="<<sum<<endl;
		return 1;
	    }

	    int pl_RR;
	    int pl_RA;
	    int pl_AA;
	    
	    if(refIsAltfromSNPfile){//the B allele is reference
		pl_RR = phredcapped( g_BB/sum );
		pl_RA = phredcapped( g_AB/sum );
		pl_AA = phredcapped( g_AA/sum );		
	    }else{//the A allele is reference
		pl_RR = phredcapped( g_AA/sum );
		pl_RA = phredcapped( g_AB/sum );
		pl_AA = phredcapped( g_BB/sum );		
	    }
	    //cerr<<pl_RR<<" "<<pl_RA<<" "<<pl_AA<<endl;
	    int minPL = MIN3( pl_RR , pl_RA , pl_AA );
	    //scaling down to 1
	    pl_RR-=minPL;
	    pl_RA-=minPL;
	    pl_AA-=minPL;
	    //cerr<<pl_RR<<" "<<pl_RA<<" "<<pl_AA<<endl;
	    pl_RR=MIN2(255,pl_RR);
	    pl_RA=MIN2(255,pl_RA);
	    pl_AA=MIN2(255,pl_AA);
	    //cerr<<pl_RR<<" "<<pl_RA<<" "<<pl_AA<<endl;
	    //toremove
	    minPL = MIN3( pl_RR , pl_RA , pl_AA );
	    if( minPL!=0){
		cerr<<"BEAGLE2GLF ERROR: the genotype likelihood for individual "<<(i/3)<<" do not sum up to 0 at "<<lineB<<"\t sum="<<sum<<endl;
		return 1;
	    }


	    SingleGL sample (uint8_t(pl_RR),
			     uint8_t(pl_RA),
			     uint8_t(pl_AA),			     
			     false);//no CpG provided
                

	    // SingleAllele sample (refCount, altCount, 0);//no CpG provided

	    arToWrite.vectorGLs->push_back(sample);
	}//for each char in lineG
	
	if(!gw->writeAlleleRecord(&arToWrite)){
	    cerr<<"beagle2glf: error writing record "<<endl;
	    exit(1);
	}	    
		
	writtenRec++;
		
	
    nextline:
	//hasDataG = (bool)getline (myfileG,lineG);
        hasDataB = (bool)getline (myfileB,lineB);
        hasDataP = (bool)getline (myfileP,lineP);
    }//while hasData


    if(rtEPO != NULL)
	delete(rtEPO);
    
    delete(gw);

    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRec<<" records, wrote "<<writtenRec<<" terminated gracefully"<<endl;

    return 0;
}

