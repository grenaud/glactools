/*
 * GlacDown
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacDown.h"

GlacDown::GlacDown(){

}

GlacDown::~GlacDown(){

}


string GlacDown::usage() const{

    
    return string("") +"glactools down [options] <acf file> [populations to downsample]\n\n"+
	"This program will downsample the allele count of the populations specified in the list.\n"+
        "Please note that it will set the alternative allele to 'N' if no population has the alternative allele and will print to STDOUT\n"+
	"If no populations are specified, it will downsample all of them\n"+
	"It will NOT downsample the root/anc unless explicitly specified in the [populations to downsample]\n"+       
	"It can be used in conjuction with bam2acf to downsample or to transform a population into a single individual using \"-c 2\"\n"+
	"\n"+
	"ex:  glactools down data.acf.gz \"Papuan,Austalian\""+"\n"+
	"will downsample the allele count for  the \"Papuan\" and \"Austalian\" individuals\n"+
	"Options:"+"\n"+	
	"\t"+"-c [count]" + "\t\t"+"Downsample exactly [count] allele (default: "+stringify(count)+")\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
	"\n"
	;
}



int GlacDown::run(int argc, char *argv[]){

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

        if(string(argv[i]) == "-c"){
            count=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }


        cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }


    if( (lastOpt != (argc-1))
	&&
	(lastOpt != (argc-2))
	){
        cerr<<"The last arguments are the  <acf file> \"popToKeep1,popToKeep2,....\"  "<<endl;
        return 1;               
    }


    string fileglac    = string(argv[lastOpt  ]);
    bool specifiedPops = (lastOpt == (argc-2));

    GlacParser gp (fileglac);
    if(gp.isGLFormat()){
        cerr<<"glactools down only works with acf files"<<endl;
        return 1;               
    }

    vector<string> g1v               ;
    vector<bool> flagsPopToAdd       ;
    string pop2keep;

    if(specifiedPops){
	pop2keep                  = string(argv[lastOpt+1]);
	g1v                       = allTokens(pop2keep,',');
	flagsPopToAdd             = vector<bool>(gp.getPopulationsNames()->size()+2,false);
        
	for(unsigned k=0;k<g1v.size();k++){
	    bool found=false;
	    for(unsigned i=0;i<(gp.getPopulationsNames()->size());i++){
		if(gp.getPopulationsNames()->at(i) == g1v[k]){
		    found=true;
		    flagsPopToAdd[i] = true;
		    break;
		} 
	    }
		
	    if(!found){
		cerr<<"Cannot find population "<<g1v[k]<<endl;
		return 1;
	    }
	}
    }else{//all except root/anc
	g1v                       = *(gp.getPopulationsNames());
	pop2keep                  = vectorToString(g1v,",");
	flagsPopToAdd             = vector<bool>(gp.getPopulationsNames()->size()+2,true);//add all	
    }



    

    uint32_t newsizepop= uint32_t( gp.getPopulationsNames()->size()-2 );
    
    GlacWriter * gw = new GlacWriter(newsizepop,
				     gp.isGLFormat(),
				     gp.isGLFormat()?1:2,
				     1,//compression threads
				     uncompressed);


    string defline =   gp.getDefline();





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
    header<<"#GLACDOWN: "<<pop2keep<<endl;

    header<<"#DOWNFILE#"<<(1)<<endl;
    header<<""<<gp.getHeaderNoSQNoDefline("#\t")<<endl;

    header<<gp.getHeaderSQ("")<<endl;
    header<<defline<<endl;
    //cout<<"test"<<header.str()<<endl;
    
    if(!gw->writeHeader(header.str())){
	cerr<<"GlacDown: error writing header "<<endl;
	exit(1);
    }
	


    AlleleRecords * arr;

    while(gp.hasData()){//TODO optimize for large pops, we do not need to read a large # of pops to extract a few

	arr = gp.getData();

	AlleleRecords  arw;
	arw.copyCoreFields(*arr);
	arw.sizePops=newsizepop;

	
	bool someoneHasAlt=false;//flag to check if someone has the alternative, otherwise, we will set the alt to N

	if(specifiedPops){
	    for(unsigned int j=0;j<arr->vectorAlleles->size();j++){
		if(flagsPopToAdd[j]){
		    arr->vectorAlleles->at(j).downsample(count);
		}

		arw.vectorAlleles->push_back(   arr->vectorAlleles->at(j));
		someoneHasAlt=someoneHasAlt || (arr->vectorAlleles->at(j).hasAlt());		
	    }
	}else{
	    for(unsigned int j=0;j<arr->vectorAlleles->size();j++){

		if(j>1){//only for all pops except the root/anc
		    arr->vectorAlleles->at(j).downsample(count);
		}
		
		arw.vectorAlleles->push_back(   arr->vectorAlleles->at(j) );
		someoneHasAlt=someoneHasAlt || (arr->vectorAlleles->at(j).hasAlt());		
	    }
	}

	if(!someoneHasAlt)
	    arw.alt='N';


	if(!gw->writeAlleleRecord(&arw)){
     	    cerr<<"GlacDown: error record "<<arw<<endl;
     	    exit(1);
     	}
    }


    delete(gw);



    return 0;
}

