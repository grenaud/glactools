/*
 * GLF2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GLF2ACF.h"

//#define DEBUGPOS 9484430

GLF2ACF::GLF2ACF(){

}

GLF2ACF::~GLF2ACF(){

}


string GLF2ACF::usage() const{

    
    return string(string("glactools") +" glf2acf [options] <glf file>  "+"\n"+
                  "\nThis program converts a GLF file containing genotype likelihoods into ACF with allele counts (prints to the stdout)\n"+    
		  
                  "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
                  
                  
                  "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"
		  );
}



int GLF2ACF::run(int argc, char *argv[]){

    int lastOpt=1;

    for(int i=1;i<(argc);i++){ 
        //cerr<<i<<"\t"<<string(argv[i])<<endl;
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
	if( string(argv[i]) == "--minPL"  ){
            minPLdiffind=destringify<int>(argv[i+1]);
            //            specifiedPL  =true;
            i++;
            continue;
        }
        cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }
    //cerr<<lastOpt<<"\t"<<(argc-1)<<" "<<minPLdiffind<<endl;
    // exit(1);
    if(lastOpt != (argc-1)){
        cerr<<"The last argument is the  <glf file>  "<<endl;
        return 1;               
    }
    string fileglf = string(argv[lastOpt]);


    GlacParser gp (fileglf);
    AlleleRecords * arr;
    AlleleRecords * arw;
    
    if(!gp.isGLFormat()){
        cerr<<"The input file must be in GLF"<<endl;
        return 1;               
    }
    
    GlacWriter * gw = new GlacWriter(gp.getSizePops(),
				     false,
				     2,
				     1,//compression threads
				     uncompressed);
    string newheader="";
    newheader+="#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    newheader+="#PG:"+programLine+"\n";;
    newheader+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";;

    newheader+="#DATE: "+getDateString()+"\n";;
    newheader+="#GLF2ACF:\n";
    newheader+=gp.getHeaderNoSQNoDefline("#\t")+"\n";
    newheader+=gp.getHeaderSQ("")+"\n";
    newheader+=gp.getDefline()+"\n";

    
    if(!gw->writeHeader(newheader)){
	cerr<<"GlacViewer: error writing header "<<endl;
	exit(1);
    }
	

        
    while(gp.hasData()){
	arr = gp.getData();
	arw = new AlleleRecords(false);
	
	arw->chr           = arr->chr;
	arw->chri          = arr->chri;
	arw->coordinate    = arr->coordinate;

	
#ifdef DEBUGPOS
	bool debugPosition=false;
	if(arw->coordinate == DEBUGPOS){
	    debugPosition=true;
	    cerr<<*arr<<endl;
	}
#endif
	//cerr<<arr->coordinate<<endl;
	arw->sizePops      = arr->sizePops;
	arw->ref           = arr->ref;
	arw->alt           = arr->alt;
	arw->vectorAlleles = new vector<SingleAllele>();
	bool foundAlt=false;
	bool hasNonZero=false;
	for(unsigned int i=0;i<arr->vectorGLs->size();i++){
	    pair<int,int> pairCount= arr->vectorGLs->at(i).returnLikelyAlleleCountForRefAlt(minPLdiffind);

#ifdef DEBUGPOS
	    if( debugPosition ){
		cerr<<i<<"\t"<<pairCount.first<<"\t"<<pairCount.second<<"\t"<<arr->vectorGLs->at(i).getIsCpg()<<endl;	    
	    }
#endif

	    hasNonZero = hasNonZero || (pairCount.first != 0) || (pairCount.second != 0);
	    foundAlt = (pairCount.second!=0) || foundAlt;
	    SingleAllele saToWrite (pairCount.first, pairCount.second, arr->vectorGLs->at(i).getIsCpg());
	    arw->vectorAlleles->push_back(saToWrite);
	}
	if(!hasNonZero)//skip positions with no information
	    continue;
	
#ifdef DEBUGPOS
	    if( debugPosition ){
		return 1;
	    }
#endif

	if(!foundAlt)
	    arw->alt           = 'N';
	
	//cout<<"write"<<*arw<<endl;
	if(!gw->writeAlleleRecord(arw)){
	    cerr<<"GlacViewer: error record "<<*arw<<endl;
	    exit(1);
	}
	//cout<<"delete"<<endl;
	delete(arw);
    }

    delete(gw);

    return 0;
}

