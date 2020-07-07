
#include "ACF2BETASCAN.h"


using namespace std;


ACF2BETASCAN::ACF2BETASCAN(){

}

ACF2BETASCAN::~ACF2BETASCAN(){

}

string ACF2BETASCAN::usage() const{
    string usage=string("glactools")+" acf2betascan  [options] <ACF file>"+
	"\nThis program will print the ACF but in betascan format\n"+
	"\nhttps://github.com/ksiewert/BetaScan\n"+

        //"The left column is the reference count. The order can be changed using --useanc and --useroot\n"+
        "Options:\n"+
	// "\t\t"+"--splitpop\t\t\tSplit pop.\n"+
	// "\t\t"+"--freq\t\t\t\tOutput frequencies\n"+
	// "\t\t"+"--onlysegsite\t\t\tUse only segregating sites\n"+
	"\t\t"+"--fold\t\t\tIgnore the ancestral/root allele to report the frequency, use minor[tab]total\n"+
	"\t\t"+"--useanc\t\t\tUse the ancestral allele to report the frequency\n"+
	"\t\t"+"--useroot\t\t\tUse the root allele to report the frequency\n"+
	"\n";
	
    return usage;
}

int ACF2BETASCAN::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }
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

	// if(string(argv[i]) == "--onlysegsite" ) {
	//     onlysegsite = true;
	//     continue;
	// }

	// if(string(argv[i]) == "--freq" ) {
	//     usefreq= true;
	//     continue;
	// }

	// if(string(argv[i]) == "--splitpop" ) {
	//     splitpop  = true;
	//     continue;
	// }

	if(string(argv[i]) == "--fold" ) {
	    fold  = true;
	    continue;
	}

	if(string(argv[i]) == "--useanc" ) {
	    useAnc  = true;
	    continue;
	}

	if(string(argv[i]) == "--useroot" ) {
	    useRoot = true;
	    continue;
	}


        cerr<<"Error unknown option #"<<argv[i]<<"#"<<endl;
        return 1;
    }


    if(useRoot && useAnc){
    	cerr<<"ACF2BETASCAN: Cannot use both the ancestor and the root"<<endl;
        return 1;           
    }
    
    if(!fold){
	if(!useRoot && !useAnc){
	    cerr<<"ACF2BETASCAN: Specify the ancestor or the root to polarize the variants"<<endl;
	    return 1;           
	}
    }else{
	if(useRoot || useAnc){
	    cerr<<"ACF2BETASCAN: Cannot specify the ancestor or the root to polarize the variants when folding"<<endl;
	    return 1;           
	}
    }
    
    string glacfile  = string(argv[lastOpt]);
    
    GlacParser gp (glacfile);

    if(!gp.isACFormat()){
        cerr<<"ACF2BETASCAN: Error the file "<<glacfile<<" should be in ACF format"<<endl;
        return 1;       
    }

    uint64_t totalRecords=0;
 
    AlleleRecords * arr;

    while(gp.hasData()){	
	arr = gp.getData();
	

	totalRecords++;

	//non-seg site, no point in looking at those
	//if(onlysegsite)
	if(!isResolvedDNA(arr->alt))
	    continue;

	unsigned int refCounter=0;
	unsigned int altCounter=0;
	vector< pair<unsigned int,unsigned int> > alleleC;

	bool rootIsRef = false;
	bool ancIsRef  = false;

	if(useRoot){
	    if(arr->vectorAlleles->at(0).getRefCount() == 0 &&
	       arr->vectorAlleles->at(0).getAltCount() == 0 ){
		continue;	    
	    }
	    
	    if(arr->vectorAlleles->at(0).getAltCount() != 0 ){
		if(arr->vectorAlleles->at(0).getRefCount() != 0 ){
		    char b=arr->vectorAlleles->at(0).generateRandomAlleleBias(arr->ref,arr->alt);
		    rootIsRef = (b==arr->ref);
		    //cerr<<"Cannot determine the root allele for "<<*arr<<endl;
		    //continue;
		    //return 1;      
		}else{	
		    rootIsRef=false;
		}
	    }else{
		rootIsRef=true;
	    }
    
	}

	if(useAnc){
	    if(arr->vectorAlleles->at(1).getRefCount() == 0 &&
	       arr->vectorAlleles->at(1).getAltCount() == 0 ){
		continue;
	    }
	    
	    if(arr->vectorAlleles->at(1).getAltCount() != 0 ){
		if(arr->vectorAlleles->at(1).getRefCount() != 0 ){
		    char b=arr->vectorAlleles->at(1).generateRandomAlleleBias(arr->ref,arr->alt);
		    ancIsRef = (b==arr->ref);
		    //cerr<<"Cannot determine the root allele for "<<*arr<<endl;
		    //return 1;      
		    //continue;
		}	
		ancIsRef=false;
	    }else{
		ancIsRef=true;
	    }
	}



	for(unsigned j=2;j<arr->vectorAlleles->size();j++){
	    //undefined site
	    refCounter+=arr->vectorAlleles->at(j).getRefCount();
	    altCounter+=arr->vectorAlleles->at(j).getAltCount();	
	}
	
	cout<<arr->coordinate<<"\t";
	
	int total=1.0;
	//if(usefreq){
	total= refCounter + altCounter;
    

	if(useRoot){
	    if(rootIsRef)
		cout<<altCounter<<"\t"<<total<<endl;
	    else//root is alt
		cout<<refCounter<<"\t"<<total<<endl;
	}

	if(useAnc){
	    if(ancIsRef)
		cout<<altCounter<<"\t"<<total<<endl;
	    else
		cout<<refCounter<<"\t"<<total<<endl;
	}

	if(fold){
	    if(refCounter<altCounter){
		cout<<refCounter<<"\t"<<total<<endl;
	    }else{
		cout<<altCounter<<"\t"<<total<<endl;
	    }
	}
	
    }//gp has data



	    

    

    cerr<<"Program ACF2BETASCAN  looked at  "<<totalRecords<<" records, terminated gracefully"<<endl;

	
    return 0;
}

