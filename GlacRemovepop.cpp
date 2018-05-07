/*
 * GlacRemovepop
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacRemovepop.h"

GlacRemovepop::GlacRemovepop(){

}

GlacRemovepop::~GlacRemovepop(){

}


string GlacRemovepop::usage() const{

    
    return string("") +"glactools removepop [options] <glac file> <populations to remove>\n"+
	"This program will remove the population specified in the list. Please note that it will set the alternative allele to 'N' if no population has the alternative allele and will print to STDOUT\n"+
	"\n"+
	"ex:  glactools removepop data.acf.gz \"Papuan,Austalian\""+"\n"+
	"will remove the \"Papuan\" and \"Austalian\" individuals and keep the others\n"+
	"Options:"+"\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\n"
	;
}



int GlacRemovepop::run(int argc, char *argv[]){

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

    if(lastOpt != (argc-2)){
        cerr<<"The last arguments are the  <glf|acf file> \"popToRemove1,popToRemove2,....\"  "<<endl;
        return 1;               
    }

    
    string fileglac                  = string(argv[lastOpt  ]);
    string pop2remove                = string(argv[lastOpt+1]);
    GlacParser gp (fileglac);

    vector<string> g1v_ = allTokens(pop2remove,',');
    vector<string> g1v;
    set<string> s1v;
    for(unsigned int i=0;i<g1v_.size();i++){
	s1v.insert( g1v_[i] );
    }
    g1v.assign( s1v.begin(), s1v.end() );

    // vector<unsigned int>   g1i;
    vector<bool> flagsPopToAdd = vector<bool>(gp.getPopulationsNames()->size()+2,true);
    //cerr<<vectorToString(g1v)<<endl;

    for(unsigned k=0;k<g1v.size();k++){
	bool found=false;
	for(unsigned i=0;i<(gp.getPopulationsNames()->size());i++){
	    if(gp.getPopulationsNames()->at(i) == g1v[k]){
		found=true;
		flagsPopToAdd[i] = false;
		break;
	    } 
	}
		
	if(!found){
	    cerr<<"Cannot find population "<<g1v[k]<<endl;
	    return 1;
	}
    }




    

    uint32_t newsizepop= uint32_t( int( gp.getSizePops() ) - int( g1v.size() ) );

    GlacWriter * gw = new GlacWriter(newsizepop,
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     1,//compression threads
				     uncompressed);



    vector<string> newdefline;
    //root and anc
    newdefline.push_back(gp.getPopulationsNames()->at(0));
    newdefline.push_back(gp.getPopulationsNames()->at(1));
    
    for(unsigned i=2;i<(gp.getPopulationsNames()->size());i++){
	if(flagsPopToAdd[i])
	    newdefline.push_back(gp.getPopulationsNames()->at(i));
    }
    string defline =    "#chr\tcoord\tREF,ALT\t"+vectorToString(newdefline,"\t");





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
    header<<"#GLACREMOVE: "<<pop2remove<<endl;

    header<<"#REMOVEPOPFILE#"<<(1)<<endl;
    header<<""<<gp.getHeaderNoSQNoDefline("#\t")<<endl;

    header<<gp.getHeaderSQ("")<<endl;
    header<<defline<<endl;
    //cout<<"test"<<header.str()<<endl;
    
    if(!gw->writeHeader(header.str())){
	cerr<<"GlacRemovePop: error writing header "<<endl;
	exit(1);
    }
	

    // cout<<"newsizepop "<<newsizepop<<endl;
    // return 0;
    AlleleRecords * arr;

    while(gp.hasData()){//TODO optimize for large pops, we do not need to read a large # of pops to extract a few

	arr = gp.getData();

	AlleleRecords  arw;
	arw.copyCoreFields(*arr);
	arw.sizePops=newsizepop;
	

	if( gp.isGLFormat() ){
	    arw.vectorGLs->push_back(    arr->vectorGLs->at(0));
	    arw.vectorGLs->push_back(    arr->vectorGLs->at(1));
	}else{
	    arw.vectorAlleles->push_back(arr->vectorAlleles->at(0));
	    arw.vectorAlleles->push_back(arr->vectorAlleles->at(1));
	}
	bool someoneHasAlt=false;//flag to check if someone has the alternative, otherwise, we will set the alt to N

			
	for(unsigned int j=2;j<arr->vectorAlleles->size();j++){

	    if(flagsPopToAdd[j]){
		if( gp.isGLFormat() ){
		    arw.vectorGLs->push_back(       arr->vectorGLs->at(j));		 
		    someoneHasAlt=someoneHasAlt || (arr->vectorGLs->at(j).hasAlt());
		}else{
		    arw.vectorAlleles->push_back(   arr->vectorAlleles->at(j));
		    someoneHasAlt=someoneHasAlt || (arr->vectorAlleles->at(j).hasAlt());
		}
	    }
	       
	}

	if(!someoneHasAlt)
	    arw.alt='N';

	// if( gp.isGLFormat() ){
	//     arw.vectorGLs->push_back(newSingleGLs);
	// }else{
	//     arw.vectorAlleles->push_back(newSingleAllele);
	// }

	if(!gw->writeAlleleRecord(&arw)){
     	    cerr<<"GlacRemovepop: error record "<<arw<<endl;
     	    exit(1);
     	}
    }


    delete(gw);



    return 0;
}

