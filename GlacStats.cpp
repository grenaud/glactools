
#include "GlacStats.h"


using namespace std;


GlacStats::GlacStats(){

}

GlacStats::~GlacStats(){

}

string GlacStats::usage() const{
    string usage=string("glactools")+" stats  [options] <ACF|GLF file>"+
	"\nThis program will print basis stats e.g. # of segregating sites\n"+
        "Options:\n"+
	// "\t\t"+"--splitpop\t\t\tSplit pop.\n"+
	// "\t\t"+"--freq\t\t\tOutput frequencies\n"+
	// "\t\t"+"--onlysegsite\t\tUse only segregating sites\n"+
	// "\t\t"+"--useanc\t\t\tUse the ancestral allele to report the frequency\n"+
	// "\t\t"+"--useroot\t\t\tUse the root allele to report the frequency\n";
	"\n";
    return usage;
}

int GlacStats::run(int argc, char *argv[]){

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

	// if(string(argv[i]) == "--useanc" ) {
	//     useAnc  = true;
	//     continue;
	// }

	// if(string(argv[i]) == "--useroot" ) {
	//     useRoot = true;
	//     continue;
	// }


        cerr<<"Error unknown option #"<<argv[i]<<"#"<<endl;
        return 1;
    }


    string glacfile  = string(argv[lastOpt]);
    
    GlacParser gp (glacfile);
    uint64_t totalRecords=0;
    uint64_t segSites    =0;
    uint64_t segSitesTS  =0;
    uint64_t segSitesTV  =0;

    uint64_t diffChr     =0;
    uint16_t chri=UINT16_MAX;

    AlleleRecords * arr;

    while(gp.hasData()){	
	arr = gp.getData();
	
	if(arr->chri != chri){
	    diffChr++;
	    chri = arr->chri ;
	}

	totalRecords++;

	//non-seg site, no point in looking at those
	
	if(!isResolvedDNA(arr->alt))
	    continue;

	segSites++;
	if(isPotentialTransition(arr->ref,arr->alt))
	    segSitesTS++;
	else
	    segSitesTV++;


	    

    }

    cout<<"Program GlacStats:"<<endl
	<<"Total records:        \t"<<totalRecords<<endl
	<<"Different chromosomes:\t"<<diffChr<<endl
	<<"Segregating sites:    \t"<<segSites<<" ("<<   100*double(segSites)/double(totalRecords) << "%)"<<endl
	<<"Segregating sites TS: \t"<<segSitesTS<<" ("<< 100*double(segSitesTS)/double(segSites) << "%)"<<endl
	<<"Segregating sites TV: \t"<<segSitesTV<<" ("<< 100*double(segSitesTV)/double(segSites)<< "%)"<<endl;

	
    return 0;
}

