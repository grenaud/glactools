
#include <iostream>
#include <fstream>
#include <memory>
#include <climits>

#include "utils.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "SimpleVCF.h"
#include "MultiVCFreader.h"
// #include "FilterVCF.h"


static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;

#include "VcfMulti2ACF.h"



VcfMulti2ACF::VcfMulti2ACF(){

}

VcfMulti2ACF::~VcfMulti2ACF(){

}

void VcfMulti2ACF::setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO){//,

			      //string & lineFromEPO){

    //lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
    lineLeftEPO=(rtEPO->readLineKS(  ));
    if(!lineLeftEPO){
	cerr<<"Error, missing data in the EPO file"<<endl;
	exit(1);
    }

    char *p, *q;
    int i;

    for (p = kstrtok(kstringPtrEPO->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
	q = (char*)aux.p;
	*q = 0;
	//cout<<i<<" >"<<p<<"<  #"<<*q<<"#"<<endl;
	if(i==0){//chr
	    epoChr                   = string(p);
	    continue;
	}
	if(i==1){//coord
	    epoCoord                 = (unsigned int)strtoul(p, NULL, 0);//strtoul(p);
	    continue;
	}
	if(i==2){//human ref
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



    // allel_anc   = fieldsEPO[3][0];//inferred ancestor
    // allel_chimp = fieldsEPO[4][0];//chimp;

}

string VcfMulti2ACF::usage() const{
    const string usage=string(string("") +" vcfm2acf [options] <vcf file> "+"\n"+
			      "\nThis program convert multi-sample VCF files into ACF (prints to the stdout)\n"+
			      // "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			      // "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
			      // "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
			      // "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
			      // "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
			      "\t"+"--fai [file]" + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: none)\n"+
			      "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	
			      "\t"+"--epo [EPO file]"       +"\t" +"Use file as EPO alignment to set the (default: none)\n"+   
			      "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   

			      "\t"+"--onlyGT"        +"\t\t" +"Do not use PL values for alleles, simply use genotypes (GT)      (default: "+booleanAsString(onlyGT)+")\n"+ 
			      //"\t"+"--epo [EPO alignment file]"       +"\t\t" +"Use file as EPO alignment   (default: none)\n"+ 			      
			      "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
			      // // "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 

			      // "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
			      // "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
			      // "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
			      // "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
			      // "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n"
			      ""
			      );
		
    //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+
    return usage;
}



//int main (int argc, char *argv[]) {
int VcfMulti2ACF::run(int argc, char *argv[]){

    // cout<<INT_MAX<<endl;
    //minPLdiffind=33;
    // int minGQcutoff=0;
    // int minMQcutoff=0;

    // int minCovcutoff=0;
    // int maxCovcutoff=1000;
    // // SetVCFFilters  * filtersVCF;
    // // bool allowCloseIndelProx = false;
    // // bool allowRepeatMasking  = false;
    // // bool allowSysErr         = false;
    // // bool allowall            = false;
    // // bool allowallMQ          = false;

    // bool ancAllele           = false;

    //double     minMapabilitycutoff =0;
    // epoFile   = "none";
    // epoFileB  = false;

    bool specifiedPL  = false;


    for(int i=1;i<(argc-1);i++){ 
	if(string(argv[i]) == "-u"){
	    uncompressed=true;
            continue;
        }

        if(string(argv[i]) == "--fai"){
            fastaIndex=string(argv[i+1]);
            i++;
            continue;
        }
            
        if( string(argv[i]) == "--onlyGT"  ){
	    onlyGT  = true;
	    //            specifiedPL  =true;
            continue;
	}
                   
        if(strcmp(argv[i],"--minPL") == 0 ){
            minPLdiffind=destringify<int>(argv[i+1]);
	    specifiedPL  =true;
            i++;
            continue;
	}

	if(strcmp(argv[i],"--epo") == 0 ){
            epoFile=string(argv[i+1]);
	    epoFileB=true;
            i++;
            continue;
	}


	cerr<<"Wrong option "<<argv[i]<<endl;
	return 1;
    }
			      
    if(argc < 2 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"Usage "<<usage()<<endl;
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

    string filenameMultiVCF = string(argv[argc-1]);
    MultiVCFreader vcfr   (filenameMultiVCF,5);



    string epoFileidx = epoFile+".tbi";

    //cerr<<"VCF file "<<(string(argv[argc-1]))<<endl;
    //cerr<<"Name pop "<<(string(argv[argc-2]))<<endl;
    // if(epoFileB)
    // 	cerr<<"EPO file "<<(epoFile)<<endl;

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
    ReadTabix * rtEPO =NULL ;
    //const kstring_t * kstringPtrEPO;

    //    string lineFromEPO;
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

    gw = new GlacWriter(int(vcfr.getPopulationNames().size()),    //gp.getSizePops(),
			false, //gp.isGLFormat(),
			2,//gp.isACFormat()?2:1,
			1,//compression threads
			uncompressed);

    stringstream header;
    header<<"#ACF"<<endl;    
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    header<<"#PG:"<<programLine<<endl;
    header<<"#GITVERSION: "<<returnGitHubVersion(argv[-1],"")<<endl;
    header<<"#DATE: "<<getDateString()<<endl;
    header<<"#VCFMULTI2ACF:"<<endl;

    map<string,uint16_t> chr2index;
    uint16_t     chrCurrentIndex=0;

    for(unsigned j=0;j<(chrFound.size());j++){
	header<< string("#SQ\tSN:")<<chrFound[j].name+"\tLN:"<<stringify(chrFound[j].length)<<"\n";
        chr2index[chrFound[j].name]=chrCurrentIndex++;
        if(chrCurrentIndex == 0xFFFF){
            cerr<<"Too many chromosomes for this build, more than 65535"<<endl;
            return 1;
        }
    }

    header<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<vectorToString(vcfr.getPopulationNames(),"\t")<<endl;
    //cout<<"#chr\tcoord\tREF,ALT\troot\t"<<namePop<<endl;

    
    if(!gw->writeHeader(header.str())){
	cerr<<"GlacViewer: error writing header "<<endl;
	return 1;
    }
    


    

    while(vcfr.hasData()){
	//cout<<"hasData"<<endl;
    	vector<SimpleVCF *> * toprint=vcfr.getMultipleData();

	if(toprint->at(0)->containsIndel() || //skip indels
	   (toprint->at(0)->getAltAlleleCount() > 1  )){ //skip tri-allelic
	    continue;
	}

	if(firstLine){
	    firstLine=false;

	    if(epoFileB){
		rtEPO = new ReadTabix( epoFile.c_str()  , 
				       epoFileidx.c_str()  , 
				       toprint->at(0)->getChr(), 
				       int(toprint->at(0)->getPosition()),INT_MAX ); //the destructor should be called automatically
		kstringPtrEPO = rtEPO->getKstringPtr();
		memset(&aux, 0, sizeof(ks_tokaux_t));
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO);
	    }
	}

	if(epoFileB){

	    if(!lineLeftEPO){
		cerr<<"Error, no data in the EPO file "<< toprint->at(0)->getChr() <<":"<< int(toprint->at(0)->getPosition()) <<endl;
		return 1;
	    }

	    if(epoChr != toprint->at(0)->getChr()){
		// cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<toprint->at(0)->getChr()<<endl;
		// return 1;
		rtEPO->repositionIterator(toprint->at(0)->getChr() , int(toprint->at(0)->getPosition()),INT_MAX);
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO);	

		if( epoChr != toprint->at(0)->getChr() ){
		    cerr<<"Error, the repositioning did not work, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<toprint->at(0)->getChr()<<endl;
		    return 1;
		}
	    }



	    while(epoCoord != toprint->at(0)->getPosition()){
		if(epoCoord > toprint->at(0)->getPosition()){
		    cerr<<"Error1, are all the sites in EPO there? Difference between coords VCF="<<(*toprint->at(0))<<"\tEPO="<<kstringPtrEPO->s<<endl;
		    return 1;
		}

		if( (toprint->at(0)->getPosition() - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO->repositionIterator(toprint->at(0)->getChr() , int(toprint->at(0)->getPosition()),INT_MAX);
		    //cout<<"repo "<<int(toprint->at(0)->getPosition())<<endl;
		}

		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO);	
	    }

	    if(epoCoord != toprint->at(0)->getPosition()){
		cerr<<"Error2, are all the sites in EPO there? Difference between coords VCF="<<( toprint->at(0)->getPosition() )<<"\tline="<<kstringPtrEPO->s<<endl;
		return 1;
	    }

	}
	
	//pair<int,int> pairCount= toprint->at(0)->returnLikelyAlleleCountForRefAlt(minPLdiffind);

	// if(pairCount.first != 0 || pairCount.second != 0 ){
	char alt=(toprint->at(0)->getAlt()=="."?'N':toprint->at(0)->getAlt()[0]);
	//string chimpString;
	//string ancString;
	SingleAllele root;
	SingleAllele anc;
	// cout<<lineFromEPO<<endl;
	// cout<<allel_chimp<<"\t"<<allel_anc<<"\t"<<epoFileB<<"\t"<<toprint->at(0)->getRef()[0]<<"\t"<<alt<<endl;
		
	if(!epoFileB){ //no epo file
	    // chimpString="0,0:0";
	    // ancString="0,0:0";
	    root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
	    anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

	    // cout<<allel_chimp<<"\t"<<toprint->at(0)->getRef()[0]<<endl;
	}else{	

	    //unresolved ancestral allele (A,C,G,T)
	    if(!isResolvedDNA(allel_chimp)){
		//chimpString="0,0:0";					
		root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
	    }
	    //resolved ancestral allele
	    else{
		if(allel_chimp == toprint->at(0)->getRef()[0]){//no diff between chimp and reference
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
		if(allel_anc == toprint->at(0)->getRef()[0]){//no diff between ancestor and reference
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
	}
	//cout<<"arToWrite"<<endl;
	AlleleRecords arToWrite (false);
	arToWrite.chri          = chr2index[toprint->at(0)->getChr()];
	arToWrite.coordinate    = toprint->at(0)->getPosition();
	arToWrite.sizePops      = int(vcfr.getPopulationNames().size());
	arToWrite.ref           = toprint->at(0)->getRef()[0];
	arToWrite.alt           = alt;		
	arToWrite.vectorAlleles = new vector<SingleAllele>  ();

	arToWrite.vectorAlleles->push_back(root);
	arToWrite.vectorAlleles->push_back(anc);
	// cout<<toprint->at(0)->getChr()<<"\t"<< toprint->at(0)->getPosition()<<"\t"<<
	// cout<<toprint->at(0)->getChr()<<"\t"<< toprint->at(0)->getPosition()<<"\t"<<
	//     toprint->at(0)->getRef()<<","<<
	//     alt<<"\t"<<
	//     chimpString<<"\t"<<
	//     ancString<<"\t";

	for(unsigned int i=0;i<toprint->size();i++){
	    pair<int,int> pairCount;
	    if( !onlyGT &&
		(toprint->at(i)->getObservedGL() ||toprint->at(i)->getObservedPL() )//has either GL or PL
	    ){
		pairCount = toprint->at(i)->returnLikelyAlleleCountForRefAlt(minPLdiffind);
	    }else{//just use GT 
		pairCount = toprint->at(i)->returnLikelyAlleleCountForRefAltJustGT();
	    }

	    SingleAllele sample (pairCount.first, pairCount.second, toprint->at(i)->isCpg() );
	    arToWrite.vectorAlleles->push_back(sample);

	    //cout<<pairCount.first<<","<<pairCount.second<<":"<<(toprint->at(i)->isCpg()?"1":"0");	

	    //cout<<toprint->at(i)->getAlleCountBasedOnGT();//	pairCount.first<<","<<pairCount.second<<":"<<(toprint->at(0)->isCpg()?"1":"0")<<endl;	
	    //if(i != (toprint->size()-1))
	    //cout<<"\t";
	}
	
	if(!gw->writeAlleleRecord(&arToWrite)){
	    cerr<<"Vcf2MultiACF: error writing header "<<endl;
	    return 1;
	}

	//	}
	//cout<<endl;
	//	}

    }

    //epoFileFP.close();
    //delete(filtersVCF);
    delete(gw);
    if(rtEPO!=NULL)
	delete(rtEPO);
    cerr<<"Program glactools vcfm2acf terminated gracefully"<<endl;
    return 0;
}

