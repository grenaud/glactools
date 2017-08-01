
#include "GlacIDXSTATS.h"


using namespace std;


GlacIDXSTATS::GlacIDXSTATS(){

}

GlacIDXSTATS::~GlacIDXSTATS(){

}

string GlacIDXSTATS::usage() const{
    string usage=string("glactools")+" idxstats  <ACF or GLF file>"+
	"\nThis program will print a tally of the number of records per chromosome\n\n";
    return usage;
}

int GlacIDXSTATS::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    string glacfile  = string(argv[argc-1]);
    string bgzf_fileidx = glacfile+".bai";
    if(!isFile(bgzf_fileidx)){
	cerr << " idxstats the index file  "<<bgzf_fileidx<<" does not seem to be present"<<endl;
	return 1;       	
    }


    
    GlacParser gp (glacfile);

    htsFile *fp = hts_open(glacfile.c_str(),"r");

    hts_idx_t *idx = sam_index_load(fp, bgzf_fileidx.c_str()); // load index
    if (idx == 0) { // index is unavailable
    	cerr<<"Cannot load index "<<bgzf_fileidx<<endl;
    	exit(1);
    }else{
    	//cerr<<"index loaded succesfully\n"; //need to have bai index
    }

    cout<<"#chromosome"<<"\t"<<"records"<<endl;	
    for(int i=0;i<gp.getNumberOfChromosomes();i++){	
	uint64_t u, v;
        hts_idx_get_stat(idx, i, &u, &v);
	cout<<gp.getChromosomeName(i)<<"\t"<<u<<endl;	
    }

	
    return 0;
}

