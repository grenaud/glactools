
#include "GlacSegsite.h"


using namespace std;


GlacSegsite::GlacSegsite(){

}

GlacSegsite::~GlacSegsite(){

}

string GlacSegsite::usage() const{
    string usage=string("glactools")+" segsite [options] <ACF file>"+
	"\nThis program will retain sites where the allele count is greater than 0 for either the reference or alternative for at least one individual\n"+
	"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	 "\t\t-ts\tKeep transitions only\n"+
	"\t\t-tv\tKeep transversions only\n";
    return usage;
}

int GlacSegsite::run(int argc, char *argv[]){

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

        if(string(argv[i]) == "-u"){
            uncompressed=true;
            continue;
        }

	if(string(argv[i]) == "-ts"){
	    onlyTS=true;
	    continue;
	}
	
	if(string(argv[i]) == "-tv"){
	    onlyTV=true;
	    continue;
	}
	
	cerr<<"Error unknown option #"<<argv[i]<<"#"<<endl;
        return 1;
    }

    if(onlyTS == true &&
       onlyTV == true ){
	cerr<<"Cannot select on transitions and transversions at the same time"<<endl;
	return 1;
    }
	       

    if(lastOpt != (argc-1)){
	cerr<<"GlacSegsite: The last argument is the <ACF file> "<<endl;
	return 1;		
    }


    string fileglf = string(argv[lastOpt]);


    GlacParser gp (fileglf);
    AlleleRecords * arr;
    if(gp.isGLFormat()){
	cerr<<"GlacSegsite: This function is only defined for  <ACF file> "<<endl;
	return 1;			
    }

    GlacWriter * gw = new GlacWriter(gp.getSizePops(),
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     1,//compression threads
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
    newheader<<"#SEGSITE:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacSegsite: error writing header "<<endl;
	exit(1);
    }
	

    uint64_t totalRecords=0;
    uint64_t keptRecords=0;
    
    
    while(gp.hasData()){
	arr = gp.getData();
	totalRecords++;


	bool haveSeenRef=false;
	bool haveSeenAlt=false;

	for(unsigned j=0;j<arr->vectorAlleles->size();j++){
	    if( (arr->vectorAlleles->at(j).getRefCount() != 0) ){
		haveSeenRef=true;
	    }
	    if( (arr->vectorAlleles->at(j).getAltCount() != 0) ){
		haveSeenAlt=true;
	    }
	}

	//print data row if have seen at least once both
	if(haveSeenRef && haveSeenAlt){
	    bool isTrans = isPotentialTransition(arr->ref,arr->alt);
	    if(onlyTV || onlyTS){

		if(onlyTV && !isTrans){		    
		    if(!gw->writeAlleleRecord(arr)){
			cerr<<"GlacSegsite: error record "<<*arr<<endl;
			exit(1);
		    }
		    keptRecords++;

		}
			
		if(onlyTS && isTrans){		    
		    if(!gw->writeAlleleRecord(arr)){
			cerr<<"GlacSegsite: error record "<<*arr<<endl;
			exit(1);
		    }
		    keptRecords++;

		}


	    }else{		    
		if(!gw->writeAlleleRecord(arr)){
		    cerr<<"GlacSegsite: error record "<<*arr<<endl;
		    exit(1);
		}
		keptRecords++;

	    }

	}

    }


    delete(gw);


    cerr<<"Program segsite wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

