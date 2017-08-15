
#include "GlacRename.h"


using namespace std;


GlacRename::GlacRename(){

}

GlacRename::~GlacRename(){

}

string GlacRename::usage() const{
    string usage=string("glactools")+" rename [options] <ACF file> \"oldpopname1,oldpopname2,...\" \"newpopname1,newpopname2,...\"  "+
	"\nThis program will rename populations  \"oldpopname1,oldpopname2,...\" with  \"newpopname1,newpopname2,...\" \n"+
	"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\n";
	//  "\t\t-ts\tKeep transitions only\n"+
	// "\t\t-tv\tKeep transversions only\n";
    return usage;
}

int GlacRename::run(int argc, char *argv[]){

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

	
	cerr<<"Error unknown option #"<<argv[i]<<"#"<<endl;
        return 1;
    }
	       

    if(lastOpt != (argc-3)){
	cerr<<"GlacRename: The last argument is the <ACF file> "<<endl;
	return 1;		
    }


    string fileglf  = string(argv[lastOpt+0]);
    string oldnames = string(argv[lastOpt+1]);
    vector<string> oldv  = allTokens(oldnames,',');
    vector<bool>   oldvf ( oldv.size(),false);

    for(unsigned int j=0;j<oldv.size();j++){
	if(oldv[j] == "root" || oldv[j] == "anc" ){
	    cerr<<"GlacRename: Cannot rename from root or anc"<<endl;
	    return 1;		
	}
    }
    string newnames = string(argv[lastOpt+2]);
    vector<string> newv = allTokens(newnames,',');
    for(unsigned int j=0;j<newv.size();j++){
	if(newv[j] == "root" || newv[j] == "anc" ){
	    cerr<<"GlacRename: Cannot rename to root or anc"<<endl;
	    return 1;		
	}
    }


    if(oldv.size() != newv.size() ){
	cerr<<"GlacRename: The size of the old names "<<oldnames<<" is not the same as the new ones "<<newnames<<endl;
	return 1;		
    }

    GlacParser gp (fileglf);
    AlleleRecords * arr;

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
    newheader<<"#RENAME:"<<oldnames<<" "<<newheader<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    string         newDefline  = gp.getDefline();
    vector<string> newDeflinev = allTokens(newDefline,'\t');    
    for(unsigned int i=5;i<newDeflinev.size();i++){
	//cerr<<newDeflinev[i]<<endl;
	for(unsigned int j=0;j<oldv.size();j++){
	    //cerr<<oldv[j]<<endl;
	    if(newDeflinev[i] == oldv[j]){		
		//cerr<<"found"<<endl;
		newDeflinev[i] = newv[j];
		oldvf[j]=true;
	    }
	}
	    
    }
    for(unsigned int j=0;j<oldvf.size();j++){
	if(!oldvf[j]){
	    cerr<<"GlacRename: pop  "<<oldv[j]<<" was not found"<<endl;
	    return 1;
	}
    }
    newheader<< vectorToString(newDeflinev,"\t") <<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacRename: error writing header "<<endl;
	return 1;
    }
	

    uint64_t totalRecords=0;
    
    
    while(gp.hasData()){
	arr = gp.getData();
	totalRecords++;


	if(!gw->writeAlleleRecord(arr)){
	    cerr<<"GlacRename: error record "<<*arr<<endl;
	    exit(1);
	}

    }


    delete(gw);


    cerr<<"Program segsite wrote  "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

