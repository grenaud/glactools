
#include "GlacUndef.h"


using namespace std;


GlacUndef::GlacUndef(){

}

GlacUndef::~GlacUndef(){

}

string GlacUndef::usage() const{
    string usage=string("glactools")+" noundef  <ACF or GLF file>"+
	"\nThis program will filter out any site where:\n"+
"\tACF: the allele count is null 0,0 for both reference and alternative\n"+
"\tGLF: the genotype likelihood does not seem defined\n"+
"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\t"+"--allowrootun\t\t\tAllow the root and ancestral to be undefined\n\n";
    return usage;
}

int GlacUndef::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }
    int lastOpt=1;
    unsigned int firstPopInd=0;
		      

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

	if(string(argv[i]) == "--allowrootun"){
	    firstPopInd=2;	    
            continue;
	}

	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;
    }


    if(lastOpt != (argc-1)){
	cerr<<"The last argument is the <ACF/GLF file> "<<endl;
	return 1;		
    }


    string fileglf = string(argv[lastOpt]);


    GlacParser gp (fileglf);
    AlleleRecords * arr;
    

    GlacWriter * gw = new GlacWriter(gp.getSizePops(),
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     uncompressed);
    stringstream newheader;
    if(gp.isGLFormat())
	newheader<<"#GLF\n";
    else
	newheader<<"#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    newheader<<"#PG:"<<programLine<<"\n";;
    newheader<<"#GITVERSION: "<<returnGitHubVersion(argv[-1],"")<<"\n";;

    newheader<<"#DATE: "<<getDateString()<<"\n";;
    newheader<<"#UNDEF:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacUndef: error writing header "<<endl;
	exit(1);
    }
	

    uint64_t totalRecords=0;
    uint64_t keptRecords=0;
    
    
    while(gp.hasData()){
	arr = gp.getData();
	totalRecords++;


	if(gp.isGLFormat()){
	    for(unsigned j=firstPopInd;j<arr->vectorGLs->size();j++){
		if( arr->vectorGLs->at(j).alleleCountNull() )
		    goto nextiterationnoundef;		
	    }
	}else{
	    for(unsigned j=firstPopInd;j<arr->vectorAlleles->size();j++){
		if( arr->vectorAlleles->at(j).alleleCountNull() )
		    goto nextiterationnoundef;
	    }

	}

	if(!gw->writeAlleleRecord(arr)){
	    cerr<<"GlacUndef: error record "<<*arr<<endl;
	    exit(1);
	}
	keptRecords++;

    nextiterationnoundef:	    
	continue;
	

    }


    delete(gw);


    cerr<<"Program undef wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

