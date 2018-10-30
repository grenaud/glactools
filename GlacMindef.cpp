
#include "GlacMindef.h"


using namespace std;


GlacMindef::GlacMindef(){

}

GlacMindef::~GlacMindef(){

}

string GlacMindef::usage() const{
    string usage=string("glactools")+" mindef [options] <ACF file>"+
	"\nThis program will filter out any site where\n"+
"\tan individual does not at least have a certain number of defined bases\n"+

"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\t"+"-min" + "\t[count]\t\t"+"Filter sites where 1 individual has less than [count] defined bases (default: "+booleanAsString(minBaseC)+")\n"+

	"\t"+"--allowrootun\t\t\tAllow the root and ancestral to be less than [count] alleles\n\n";
    return usage;
}

int GlacMindef::run(int argc, char *argv[]){

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

	if(string(argv[i]) == "-min"){
	    minBaseC=destringify<int>(argv[i+1]);	    
	    i++;
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
    if(gp.isGLFormat()){
	cerr<<"The input file must be in GLF"<<endl;
	return 1;               
    }
    AlleleRecords * arr;
    

    GlacWriter * gw = new GlacWriter(gp.getSizePops(),
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     1,//compression threads
				     uncompressed);
    stringstream newheader;
    newheader<<"#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    newheader<<"#PG:"<<programLine<<"\n";;
    newheader<<"#GITVERSION: "<<returnGitHubVersion(argv[-1],"")<<"\n";;

    newheader<<"#DATE: "<<getDateString()<<"\n";;
    newheader<<"#MINDEF:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacMindef: error writing header "<<endl;
	exit(1);
    }
	

    uint64_t totalRecords=0;
    uint64_t keptRecords=0;
    
    
    while(gp.hasData()){
	arr = gp.getData();
	totalRecords++;


	for(unsigned j=firstPopInd;j<arr->vectorAlleles->size();j++){
	    if( arr->vectorAlleles->at(j).getTotalCount() < minBaseC )
		goto nextiterationnomindef;
	}

	if(!gw->writeAlleleRecord(arr)){
	    cerr<<"GlacMindef: error writing record "<<*arr<<endl;
	    exit(1);
	}
	keptRecords++;

    nextiterationnomindef:	    
	continue;
	

    }


    delete(gw);


    cerr<<"Program mindef wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

