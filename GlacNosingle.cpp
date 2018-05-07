
#include "GlacNosingle.h"


using namespace std;


GlacNosingle::GlacNosingle(){

}

GlacNosingle::~GlacNosingle(){

}

string GlacNosingle::usage() const{
    string usage=string("glactools")+" sharing [options] <ACF file> "+
	"\nThis will only retain sites where 2 individuals or more have a variant.\n"+
	"which could be the reference or alternative base\n"
	"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\t"+"-r" + "\t\t\t"+"Include root/anc as potential individuals (default: "+booleanAsString(includeRootAnc)+")\n";

    return usage;
}

int GlacNosingle::run(int argc, char *argv[]){

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
            uncompressed   = true;
            continue;
        }

        if(string(argv[i]) == "-r"){
            includeRootAnc = true;
            continue;
        }

	cerr<<"Error unknown option #"<<argv[i]<<"#"<<endl;
        return 1;
    }
	       

    if(lastOpt != (argc-1)){
	cerr<<"GlacNosingle: The last arguments are the <ACF file> "<<endl;
	return 1;		
    }


    string fileacf = string(argv[lastOpt]);

	    


    GlacParser gp (fileacf);
    AlleleRecords * arr;
    if(gp.isGLFormat()){
	cerr<<"GlacNosingle: This function is only defined for  <ACF file> "<<endl;
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
    newheader<<"#NOSINGLE:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacNosingle: error writing header "<<endl;
	exit(1);
    }
	

    uint64_t totalRecords=0;
    uint64_t keptRecords=0;
    
    unsigned jFirst;
    if(includeRootAnc){
	jFirst = 0;
    }else{
	jFirst = 2;
    }
    while(gp.hasData()){
	arr = gp.getData();

	totalRecords++;
	
	if(arr->ref == 'N' ||
	   arr->alt == 'N' )	       
	    continue;


	int refCountg=0;
	int altCountg=0;

	

	for(unsigned j=jFirst;j<arr->vectorAlleles->size();j++){
		    

	    if(arr->vectorAlleles->at(j).getRefCount() == 0 &&
	       arr->vectorAlleles->at(j).getAltCount() == 0 ){
		continue;
	    }

	    if(arr->vectorAlleles->at(j).getRefCount() != 0 &&
	       arr->vectorAlleles->at(j).getAltCount() == 0 ){
		refCountg++;
		continue;
	    }
	    
	    if(arr->vectorAlleles->at(j).getRefCount() == 0 &&
	       arr->vectorAlleles->at(j).getAltCount() != 0 ){
		altCountg++;
		continue;
	    }
	    
	}
	//cerr<<arr->coordinate<<" "<<refCountg<<" "<<altCountg<<endl;		
	if (refCountg <= 1 ){ //only one has the ref
	    continue;
	}

	if (altCountg <= 1 ){ //only one has the alt
	    continue;
	}

	//print data

	if(!gw->writeAlleleRecord(arr)){
	    cerr<<"GlacNosingle: error record "<<*arr<<endl;
	    exit(1);
	}
	keptRecords++;
		

    }


    delete(gw);


    cerr<<"Program nosingle wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

