
#include "ACF2EIGENSTRAT.h"


using namespace std;


ACF2EIGENSTRAT::ACF2EIGENSTRAT(){

}

ACF2EIGENSTRAT::~ACF2EIGENSTRAT(){

}

string ACF2EIGENSTRAT::usage() const{
    string usage=string("glactools")+" acf2eigenstrat  [options] <ACF file> [out prefix]"+
	"\nThis program takes a ACF file and prints the output as EIGENSTRAT\n\n"+
	"\tOptions\n"+			
        "\t\t"+"--withanc"+"\t"+"Print the anc  (Default: "+boolStringify(printAnc)+" )\n"+
        "\t\t"+"--withroot"+"\t"+"Print the root (Default: "+boolStringify(printRoot)+" )\n"+
	"";

    return usage;
}

int ACF2EIGENSTRAT::run(int argc, char *argv[]){



    if(argc < 2  ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    //all but last 2
    int lastOpt=1;

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

	if( string(argv[i]) == "--withroot"){
	    printRoot=true;
	    continue;
	}

	if( string(argv[i]) == "--withanc"){
	    printAnc =true;
	    continue;
	}


	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    string glacfile  = string(argv[lastOpt]);
    
    string genoFile   = string(argv[lastOpt+1])+".geno";
    string snpFile    = string(argv[lastOpt+1])+".snp";
    string indFile    = string(argv[lastOpt+1])+".ind";

    GlacParser gp (glacfile);
    if(!gp.isACFormat()){
	cerr<<"ACF2EIGENSTRAT: Error the file "<<glacfile<<" should be in ACF format"<<endl;
        return 1;	
    }

    ofstream genoFileS;
    ofstream snpFileS; 
    ofstream indFileS; 
	
    genoFileS.open(genoFile.c_str(), ios::out);
    snpFileS.open(snpFile.c_str(),   ios::out);
    indFileS.open(indFile.c_str(),   ios::out);

    
    if (!genoFileS.good()){       cerr << "Unable to open file "<<genoFile<<endl;       return 1;     }
    if (!snpFileS.good()){        cerr << "Unable to open file "<<snpFile<<endl;        return 1;     }
    if (!indFileS.good()){        cerr << "Unable to open file "<<indFile<<endl;        return 1;     }
    
    

    AlleleRecords * record;
 


    unsigned int firstIndex=2;
    if(printRoot)
	indFileS<<gp.getPopulationsNames()->at(0)<<"\tU\t"<<gp.getPopulationsNames()->at(0)<<endl;
    if(printAnc)
	indFileS<<gp.getPopulationsNames()->at(1)<<"\tU\t"<<gp.getPopulationsNames()->at(1)<<endl;

    for(unsigned int i=firstIndex;i<gp.getPopulationsNames()->size();i++){
	indFileS<<gp.getPopulationsNames()->at(i)<<"\tU\t"<<gp.getPopulationsNames()->at(i)<<endl;
	
    } 
    indFileS.close();



    unsigned int counter=0;
    while(gp.hasData()){
	record = gp.getData();
	
	if(!isResolvedDNA(record->alt))
	    continue;

	snpFileS<<"snp#"<<(counter++)<<"\t"<<record->chr<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	unsigned int firstIndex=2;
	if(printRoot)
	    genoFileS<<record->vectorAlleles->at(0).printEIGENSTRAT();	   
	if(printAnc)
	    genoFileS<<record->vectorAlleles->at(1).printEIGENSTRAT();	   

	for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
	    genoFileS<<record->vectorAlleles->at(i).printEIGENSTRAT();	   
	} 
	genoFileS<<endl;

    }
    genoFileS.close();
    snpFileS.close();
        
    cerr<<"Program "<<argv[0]<<" terminated gracefully, wrote "<<counter<<" sites"<<endl;
	
    return 0;
}

