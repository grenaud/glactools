/*
 * GlacMeld
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacMeld.h"

GlacMeld::GlacMeld(){

}

GlacMeld::~GlacMeld(){

}


string GlacMeld::usage() const{

    
    return string("") +"glactools glacmeld [options] <glac file> \"popToMerge1,popToMerge2,....\" \"newid\"\n"+
	"This program will merge different specified populations into a single one and will print to STDOUT\n"+
	"\n"+
	"ex:  glactools glacmeld data.acf.gz \"Papuan,Austalian\" \"oceanians\""+"\n"+
	"\n"+
	"Options:"+"\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+

	"\t--keep\t\t\tKeep the original populations in the output (Default "+boolStringify(keepOrig)+" )\n";
    ;
}



int GlacMeld::run(int argc, char *argv[]){

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

        if(string(argv[i]) == "-k"){
            keepOrig=true;
            continue;
        }

        cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    if(lastOpt != (argc-3)){
        cerr<<"The last arguments are the  <acf file> \"popToMerge1,popToMerge2,....\" \"newid\" "<<endl;
        return 1;               
    }

    
    string fileglac                  = string(argv[lastOpt  ]);
    string pop2mergeString           = string(argv[lastOpt+1]);
    string mergedpopName             = string(argv[lastOpt+2]); 

    vector<string> tomerge = allTokens( pop2mergeString,',');
    set<unsigned int> indexPopTomerge;

    GlacParser gp (fileglac);
    if(	gp.isGLFormat()){
        cerr<<"meld cannot merge genotype likelihood files GLF"<<endl;
        return 1;               
    }

    

    uint32_t newsizepop= 0;
    if(keepOrig){ 
	newsizepop = gp.getSizePops()+1 ; 
    }else{ 
	newsizepop = uint32_t(int(gp.getSizePops())-int(tomerge.size())+1) ; 
    }
    GlacWriter * gw = new GlacWriter(newsizepop,
				     false,
				     2,
				     1,//compression threads
				     uncompressed);

    for(unsigned int i=0;i<gp.getPopulationsNames()->size();i++){
	for(unsigned int j=0;j<tomerge.size();j++){
	    if( tomerge[j] == gp.getPopulationsNames()->at(i) ){
		indexPopTomerge.insert(i);
	    }
	}
    }

    if(indexPopTomerge.size() != tomerge.size()){
	cerr<<"Error: some of the populations to merge "<<pop2mergeString<<" were not found"<<endl;
	return 1;
    }
    string defline =    "#chr\tcoord\tREF,ALT\t";
    for(unsigned int i=0;i<gp.getPopulationsNames()->size();i++){
	if( indexPopTomerge.find(i) ==   indexPopTomerge.end() ){
	    defline+=gp.getPopulationsNames()->at(i)+"\t";
	}else{
	    if(keepOrig)
		defline+=gp.getPopulationsNames()->at(i)+"\t";
	    //skip
	}
    }

    defline+=mergedpopName+"\n";



    stringstream header;
    header<<"#ACF"<<endl;    	

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    header<<"#PG:"<<programLine<<endl;
    header<<"#GITVERSION: "<<returnGitHubVersion(argv[-1],"")<<endl;
    header<<"#DATE: "<<getDateString()<<endl;
    header<<"#GLACMELD: "<<pop2mergeString<<" "<<mergedpopName<<endl;

    header<<"#MELDFILE#"<<(1)<<endl;
    header<<""<<gp.getHeaderNoSQNoDefline("#\t")<<endl;

    header<<gp.getHeaderSQ("")<<endl;
    header<<defline<<endl;

    //cout<<header.str()<<endl;
    
    if(!gw->writeHeader(header.str())){
	cerr<<"GlacMeld: error writing header "<<endl;
	exit(1);
    }
	

    //cout<<"newsizepop "<<newsizepop<<endl;
    AlleleRecords * arr;

    while(gp.hasData()){
	//cout<<endl<<"data"<<endl;
	arr = gp.getData();
	// cout<<*arr<<endl;
	// cout<<arr->chr<<"\t"<<arr->coordinate<<"\t"<<arr->ref<<","<<arr->alt<<"\t"<<endl;

	AlleleRecords  arw;
	arw.copyCoreFields(*arr);
	arw.sizePops=newsizepop;
	SingleAllele newSingleAllele;
	// newSingleAllele.refCount=0;
	// newSingleAllele.altCount=0;
	// newSingleAllele.isCpg=false;
	
	for(unsigned int i=0;i<arr->vectorAlleles->size();i++){

	    if(indexPopTomerge.find(i) ==   indexPopTomerge.end() ){//print normally if not found
		arw.vectorAlleles->push_back(arr->vectorAlleles->at(i));
	    }else{
		if(keepOrig){
		    //cout<<test->vectorAlleles->at(i) <<"\t";
		    arw.vectorAlleles->push_back(arr->vectorAlleles->at(i));
		}
		newSingleAllele+=arr->vectorAlleles->at(i);
	    }
	}
	
	//cout<<newSingleAllele <<endl; 
	arw.vectorAlleles->push_back(newSingleAllele);

	if(!gw->writeAlleleRecord(&arw)){
     	    cerr<<"GlacMeld: error record "<<arw<<endl;
     	    exit(1);
     	}
    }

    // while(gp.hasData()){
    // 	arr = gp.getData();
    // 	arw = new AlleleRecords(false);
	
    // 	arw->chr           = arr->chr;
    // 	arw->chri          = arr->chri;
    // 	arw->coordinate    = arr->coordinate;
    // 	//cout<<arr->coordinate<<endl;
    // 	arw->sizePops      = arr->sizePops;
    // 	arw->ref           = arr->ref;
    // 	arw->alt           = arr->alt;
    // 	arw->vectorAlleles = new vector<SingleAllele>();
	
    // 	for(unsigned int i=0;i<arr->vectorGLs->size();i++){
    // 	    pair<int,int> pairCount= arr->vectorGLs->at(i).returnLikelyAlleleCountForRefAlt(minPLdiffind);
    // 	    //cout<<i<<"\t"<<pairCount.first<<"\t"<<pairCount.second<<"\t"<<arr->vectorGLs->at(i).getIsCpg()<<endl;
    // 	    SingleAllele saToWrite (pairCount.first, pairCount.second, arr->vectorGLs->at(i).getIsCpg());
    // 	    arw->vectorAlleles->push_back(saToWrite);
    // 	}
    // 	//cout<<"write"<<*arw<<endl;
    // 	if(!gw->writeAlleleRecord(arw)){
    // 	    cerr<<"GlacViewer: error record "<<*arw<<endl;
    // 	    exit(1);
    // 	}
    // 	//cout<<"delete"<<endl;
    // 	delete(arw);
    // }

    delete(gw);



    return 0;
}

