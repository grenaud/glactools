/*
 * GLF2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GLF2ACF.h"

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
	if( string(argv[i]) == "--minPL"  ){
            minPLdiffind=destringify<int>(argv[i+1]);
            //            specifiedPL  =true;
            i++;
            continue;
        }

    }

    if(lastOpt != (argc-1)){
        cerr<<"The last argument is the  <glf file>  "<<endl;
        return 1;               
    }
    string fileglf = string(argv[lastOpt]);


    GlacParser gp (fileglf);
    AlleleRecords * arr;
    AlleleRecords * arw;
    

    GlacWriter * gw = new GlacWriter(gp.getSizePops(),
				     false,
				     2,
				     uncompressed);
    string newheader="";
    newheader+="#GLF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

    newheader+="#PG:"+programLine+"\n";;
    newheader+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";;

    newheader+="#DATE: "+getDateString()+"\n";;
    newheader+="#GLF2ACF:";
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
	//cout<<arr->coordinate<<endl;
	arw->sizePops      = arr->sizePops;
	arw->ref           = arr->ref;
	arw->alt           = arr->alt;
	arw->vectorAlleles = new vector<SingleAllele>();
	
	for(unsigned int i=0;i<arr->vectorGLs->size();i++){
	    pair<int,int> pairCount= arr->vectorGLs->at(i).returnLikelyAlleleCountForRefAlt(minPLdiffind);
	    //cout<<i<<"\t"<<pairCount.first<<"\t"<<pairCount.second<<"\t"<<arr->vectorGLs->at(i).getIsCpg()<<endl;
	    SingleAllele saToWrite (pairCount.first, pairCount.second, arr->vectorGLs->at(i).getIsCpg());
	    arw->vectorAlleles->push_back(saToWrite);
	}
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

