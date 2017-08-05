
#include "Glac2FREQSPEC.h"


using namespace std;


Glac2FREQSPEC::Glac2FREQSPEC(){

}

Glac2FREQSPEC::~Glac2FREQSPEC(){

}

string Glac2FREQSPEC::usage() const{
    string usage=string("glactools")+" freqspec  [options] <ACF file>"+
	"\nThis program will print the number of observed alleles for the reference and alternative alleles\n"+
        "The left column is the reference count. The order can be changed using --useanc and --useroot\n"+
        "Options:\n"+
	"\t\t"+"--splitpop\t\t\tSplit pop.\n"+
	"\t\t"+"--freq\t\t\tOutput frequencies\n"+
	"\t\t"+"--onlysegsite\t\tUse only segregating sites\n"+
	"\t\t"+"--useanc\t\t\tUse the ancestral allele to report the frequency\n"+
	"\t\t"+"--useroot\t\t\tUse the root allele to report the frequency\n";

    return usage;
}

int Glac2FREQSPEC::run(int argc, char *argv[]){

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

	if(string(argv[i]) == "--onlysegsite" ) {
	    onlysegsite = true;
	    continue;
	}

	if(string(argv[i]) == "--freq" ) {
	    usefreq= true;
	    continue;
	}

	if(string(argv[i]) == "--splitpop" ) {
	    splitpop  = true;
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
	cerr<<"Glac2FREQSPEC: Cannot use both the ancestor and the root"<<endl;
        return 1;           
    }
    string glacfile  = string(argv[argc-1]);
    
    GlacParser gp (glacfile);
    uint64_t totalRecords=0;
 
    AlleleRecords * arr;

    while(gp.hasData()){	
	arr = gp.getData();
	

	totalRecords++;

	//non-seg site, no point in looking at those
	if(onlysegsite)
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
		    cerr<<"Cannot determine the root allele for "<<*arr<<endl;
		    return 1;      
		}	
		rootIsRef=false;
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
		    cerr<<"Cannot determine the root allele for "<<*arr<<endl;
		    return 1;      
		}	
		ancIsRef=false;
	    }else{
		ancIsRef=true;
	    }
	}



	if(splitpop){
	    for(unsigned j=2;j<arr->vectorAlleles->size();j++){
		alleleC.push_back( make_pair<unsigned int,unsigned int>( arr->vectorAlleles->at(j).getRefCount(),
									 arr->vectorAlleles->at(j).getAltCount() ) );
	    }

	    cout<<arr->coordinate;
	    bool flipAlt=false;
	    if(useRoot){
		if(!rootIsRef)
		    flipAlt=true;
	    }
	    
	     if(useAnc){
		if(!ancIsRef)
		    flipAlt=true;
	     }
	     
	     if(flipAlt){
		 
		 for(unsigned j=0;j<alleleC.size();j++){
		     double total=1.0;
		     if(usefreq)
			 total= double(alleleC[j].first) + double(alleleC[j].second);
		     cout<<"\t"<<alleleC[j].second/total<<"\t"<<alleleC[j].first/total;

		 }
		 
	     }else{
		 for(unsigned j=0;j<alleleC.size();j++){
		     double total=1.0;
		     if(usefreq)
			 total= double(alleleC[j].first) + double(alleleC[j].second);
		     cout<<"\t"<<alleleC[j].first/total<<"\t"<<alleleC[j].second/total;
		     
		 }
	     }
	     cout<<endl;
	}else{
	    for(unsigned j=2;j<arr->vectorAlleles->size();j++){
		//undefined site
		refCounter+=arr->vectorAlleles->at(j).getRefCount();
		altCounter+=arr->vectorAlleles->at(j).getAltCount();	
	    }

	    cout<<arr->chr<<"\t"<<arr->coordinate<<"\t";
	    double total=1.0;
	    if(usefreq){
		total= double(refCounter) + double(altCounter);
	    }

	    if(useRoot){
		if(rootIsRef)
		    cout<<double(refCounter)/total<<"\t"<<double(altCounter)/total<<endl;
		else
		    cout<<double(altCounter)/total<<"\t"<<double(refCounter)/total<<endl;
		continue;
	    }

	    if(useAnc){
		if(ancIsRef)
		    cout<<double(refCounter)/total<<"\t"<<double(altCounter)/total<<endl;
		else
		    cout<<double(altCounter)/total<<"\t"<<double(refCounter)/total<<endl;
		continue;
	    }

	    cout<<double(refCounter)/total<<"\t"<<double(altCounter)/total<<endl;
	}



	    

    }

    cerr<<"Program Glac2FREQSPEC  looked at  "<<totalRecords<<" records, terminated gracefully"<<endl;

	
    return 0;
}

