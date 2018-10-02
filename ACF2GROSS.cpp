
#include "ACF2GROSS.h"


using namespace std;


ACF2GROSS::ACF2GROSS(){

}

ACF2GROSS::~ACF2GROSS(){

}

string ACF2GROSS::usage() const{


    string usage=string("")+"acf2gross  [options]  <ACF file>\n"+
	"This program takes an ACF matrix and prints the segregating sites in GRoSS format (https://github.com/FerRacimo/GRoSS)\n"+
	"\n"+
	"Options:\n"+
	"\t--noroot\t\tDo not print the root/anc  (Default "+boolStringify(noroot)+" )\n"+
	"\t--justtransv\t\tOnly allow transversions (Default "+boolStringify(limitToTransversions)+" )\n"+
	"\t--noprivate\t\tDo not allow private mutations (Default "+boolStringify(noprivate)+" )\n"+
	"\t--anc\t\tPrint ancestral value, not the root (Default "+boolStringify(printAnc)+" )\n";
    	   	
    return usage;
}

int ACF2GROSS::run(int argc, char *argv[]){



    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    int lastOpt=1;    
    for(int i=1;i<(argc);i++){ 
	
	if((string(argv[i]) == "-")  ){
	    lastOpt=i;
	    break;          
	}
	
	if(string(argv[i])[0] != '-' ){
	    lastOpt=i;
	    break;
	}       	

        if( string(argv[i]) == "--noroot"  ){
	    noroot=true;
            continue;
        }

        if( string(argv[i]) == "--anc"  ){
	    printAnc=true;
            continue;
        }

        if( string(argv[i]) == "--justtransv"  ){
	    limitToTransversions=true;
            continue;
        }

        if( string(argv[i]) == "--noprivate"  ){
	    noprivate=true;
            continue;
        }

	cerr<<"Error unknown option "<<argv[i]<<endl;
	exit(1);
    }

    if(printAnc && !noroot){
	cerr<<"ACF2GROSS Error: cannot specify --noroot and --anc "<<endl;
	exit(1);
    }

    GlacParser gp (argv[lastOpt]);

    vector<string> toprintPop;
    unsigned jinit=0;
    if(noroot){
	jinit=2;
    }

    for(unsigned j=jinit;j<gp.getPopulationsNames()->size();j++){

	if(j == 0){
	    if(printAnc)
		continue;
	    toprintPop.push_back(gp.getPopulationsNames()->at(j));
	    continue;
	}
	
	if(j == 1){
	    if(!printAnc)
		continue;
	    toprintPop.push_back(gp.getPopulationsNames()->at(j));
	    continue;
	}
	toprintPop.push_back(gp.getPopulationsNames()->at(j));
    }

    cout<<"CHROM\tPOS\tSNPID\t"<<vectorToString( toprintPop,"\t")<<endl;
    AlleleRecords * test;
    unsigned int totalRecords=0;
    unsigned int keptRecords=0;

    
    while(gp.hasData()){
	//cout<<"data"<<endl;
	test = gp.getData();
	totalRecords++;
	if( (totalRecords%10000000) == 0){
	    cerr<<"Looked at "<<totalRecords<<" records @ "<<test->chr<<":"<<test->coordinate<<endl;
	}
	
	if(test->alt == 'N'){
	    continue;
	}

	
	if(limitToTransversions){
	    //skip potential transitions
	    if(isPotentialTransition(test->ref,test->alt))
		continue;
	}


	string toprint=stringify(test->chr)+"\t"+stringify(test->coordinate)+"\tSNP"+stringify(keptRecords)+"\t";
	int counterIndRef=0;
	int counterIndAlt=0;

	for(unsigned j=jinit;j<test->vectorAlleles->size();j++){

	    if(j == 0){
		if(printAnc)
		    continue;
	    }

	    if(j == 1){
		if(!printAnc)
		    continue;
	    }



	    if( (test->vectorAlleles->at(j).getRefCount() == 0) && 
		(test->vectorAlleles->at(j).getAltCount() == 0) ){
		goto nextiteration;
	    }

	    if(test->vectorAlleles->at(j).getRefCount() != 0)
		counterIndRef++;
	    if(test->vectorAlleles->at(j).getAltCount() != 0)
		counterIndAlt++;
	    toprint+=stringify(test->vectorAlleles->at(j).getRefCount())+","+stringify(test->vectorAlleles->at(j).getAltCount());

	    if(j!=  (test->vectorAlleles->size()-1)){
		toprint+="\t";
	    }	    
	}

	//cout<<endl;
	if(noprivate){
	    //each allele was observed more than once
	    if(counterIndRef >1 &&
	       counterIndAlt >1 ){
		cout<<toprint<<endl;		
		keptRecords++;
	    }else{
		continue;
	    }      
	}else{
	    cout<<toprint<<endl;
	    keptRecords++;
	}
    nextiteration:
	continue;
    }

    cerr<<"Program "<<argv[0]<<" wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

