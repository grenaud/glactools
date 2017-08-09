/*
 * GlacUsePopAsRootAnc
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacUsePopAsRootAnc.h"

GlacUsePopAsRootAnc::GlacUsePopAsRootAnc(){

}

GlacUsePopAsRootAnc::~GlacUsePopAsRootAnc(){

}


string GlacUsePopAsRootAnc::usage() const{

    
    return string("") +"glactools removepop [options] <glac file> <pop to use as root> <pop to use as anc>\n"+
	"This program will use specified populations <pop to use as root> as root and  <pop to use as anc> ancestor and produce lines with only those two populations. It will print to STDOUT\n"+
	"\n"+
	"ex:  glactools usepopsrootanc data.acf.gz  Pop1 Pop2"+"\n"+
	"will use the Pop1 and Pop2 as the root and anc\n"+
	"Options:"+"\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\n"
	;
}



int GlacUsePopAsRootAnc::run(int argc, char *argv[]){

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

        if(string(argv[i]) == "-u"){
            uncompressed=true;
            continue;
        }


        cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    if(lastOpt != (argc-3)){
        cerr<<"The last arguments are the  <glf|acf file> popToUseAsRoot popToUseAsAnc "<<endl;
        return 1;               
    }

    
    string fileglac          = string(argv[lastOpt  ]);
    string poproot           = string(argv[lastOpt+1]);
    string popanc            = string(argv[lastOpt+1]);

    GlacParser gp (fileglac);    

    unsigned int rootidx =0;
    unsigned int ancidx  =1;
    bool rootflg = false;
    bool ancflg  = false;

    for(unsigned int i=0;i<gp.getPopulationsNames()->size();i++){

        if(gp.getPopulationsNames()->at(i) ==  poproot){
            rootidx =i;
            rootflg=true;
        }

        if(gp.getPopulationsNames()->at(i) ==  popanc){
            ancidx =i;
            ancflg=true;
        }

    }

    if(!rootflg ){
        cerr<<"GlacUsePopAsRootAnc: Error: the root pop was not found"<<endl;
        return 1;
    }

    if(!ancflg){
        cerr<<"GlacUsePopAsRootAnc: Error: the root pop was not found"<<endl;
        return 1;
    }


    GlacWriter * gw = new GlacWriter(0,
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     1,//compression threads
				     uncompressed);



    stringstream header;
    if(gp.isGLFormat())
        header<<"#GLF"<<endl;           
    else
        header<<"#ACF"<<endl;           

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    header<<"#PG:"<<programLine<<endl;
    header<<"#GITVERSION: "<<returnGitHubVersion(argv[-1],"")<<endl;
    header<<"#DATE: "<<getDateString()<<endl;
    header<<"#USEPOPASROOTANC: "<<poproot<<" "<<popanc<<endl;
    header<<"#USEPOPASROOTANC#"<<(1)<<endl;
    header<<""<<gp.getHeaderNoSQNoDefline("#\t")<<endl;

    header<<gp.getHeaderSQ("")<<endl;
    header<<"#chr\tcoord\tREF,ALT\troot\tanc"<<endl;

    
    if(!gw->writeHeader(header.str())){
	cerr<<"GlacRemovePop: error writing header "<<endl;
	exit(1);
    }
	

    AlleleRecords * arr;

    while(gp.hasData()){//TODO optimize for large pops, we do not need to read a large # of pops to extract a few

	arr = gp.getData();

	AlleleRecords  arw;
	arw.copyCoreFields(*arr);
	arw.sizePops=0;
	

	if( gp.isGLFormat() ){
	    arw.vectorGLs->push_back(    arr->vectorGLs->at(rootidx));
	    arw.vectorGLs->push_back(    arr->vectorGLs->at(ancidx));
	}else{
	    arw.vectorAlleles->push_back(arr->vectorAlleles->at(rootidx));
	    arw.vectorAlleles->push_back(arr->vectorAlleles->at(ancidx));
	}

	if(!gw->writeAlleleRecord(&arw)){
     	    cerr<<"GlacUsePopAsRootAnc: error record "<<arw<<endl;
     	    return 1;
     	}
    }


    delete(gw);



    return 0;
}

