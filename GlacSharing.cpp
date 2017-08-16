
#include "GlacSharing.h"


using namespace std;


GlacSharing::GlacSharing(){

}

GlacSharing::~GlacSharing(){

}

string GlacSharing::usage() const{
    string usage=string("glactools")+" sharing [options] <ACF file>  <group 1> <group 2>"+
	"\nThis will only retain sites where every individuals in population group 1 share the same allele(s) as every individual in population group 2.\n"+
	"The individuals/populations in <group 1> and <group 2> are comma-separated\n"
	"In other words:\n"+
	"If the individuals in group 1 have the reference allele, it must be observed in group2\n"+
	"If the individuals in group 1 have the alternative allele, it must be observed in group2\n"+
	"Example: glactools sharing input.acf.gz \"CEU,TSI\" \"CHB,JPT\"\n"+
	"It requires that the allele count for every individual for both groups be non-zero.\n"+
	"A random allele is picked (biased for allele count) for heterozygous (or variable in a population)  position so do not be surprised if you get different outputs every time\n"+
	"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n";
    return usage;
}

int GlacSharing::run(int argc, char *argv[]){

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
	cerr<<"GlacSharing: The last arguments are the <ACF file> "<<endl;
	return 1;		
    }


    string fileglf = string(argv[lastOpt]);
    string g1 = string(argv[lastOpt+1]);
    string g2 = string(argv[lastOpt+2]);
    vector<string> g1v = allTokens(g1,',');
    vector<string> g2v = allTokens(g2,',');
    vector<unsigned int>   g1i;
    vector<unsigned int>   g2i;

	    


    GlacParser gp (fileglf);
    AlleleRecords * arr;
    if(gp.isGLFormat()){
	cerr<<"GlacSharing: This function is only defined for  <ACF file> "<<endl;
	return 1;			
    }


    for(unsigned k=0;k<g1v.size();k++){
	for(unsigned j=0;j<g2v.size();j++){
	    if(g1v[k] ==  g2v[j]){
		cerr<<"Cannot specify the same population ("<<g1v[k]<<") in both groups"<<endl;
		return 1;
	    }
	}
    }	    

    for(unsigned k=0;k<g1v.size();k++){
	bool found=false;
	for(unsigned i=0;i<(gp.getPopulationsNames()->size());i++){
	    if(gp.getPopulationsNames()->at(i) == g1v[k]){
		g1i.push_back(i);
		found=true;
		break;
	    } 
	}
	if(!found){
	    cerr<<"Cannot find population "<<g1v[k]<<endl;
	    return 1;
	}
    }

    for(unsigned k=0;k<g2v.size();k++){
	bool found=false;
	for(unsigned i=0;i<gp.getPopulationsNames()->size();i++){
	    if(gp.getPopulationsNames()->at(i) == g2v[k]){
		g2i.push_back(i);
		found=true;
		break;
	    } 
	}
	if(!found){
	    cerr<<"Cannot find population "<<g2v[k]<<endl;
	    return 1;
	}
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
    newheader<<"#SHARING:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacSharing: error writing header "<<endl;
	exit(1);
    }
	

    uint64_t totalRecords=0;
    uint64_t keptRecords=0;
    
    
    while(gp.hasData()){
	arr = gp.getData();
	totalRecords++;
	int refCountg1=0;
	int refCountg2=0;
	int altCountg1=0;
	int altCountg2=0;



	for(unsigned j=0;j<arr->vectorAlleles->size();j++){
		    
		    
	    for(unsigned k=0;k<g1i.size();k++){
		if(j==g1i[k]){
		    //skip undefined sites 
		    if(arr->vectorAlleles->at(j).getRefCount() == 0 && 
		       arr->vectorAlleles->at(j).getAltCount() == 0 )
			goto nextiterationshare;			       

		    if(arr->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
		       arr->vectorAlleles->at(j).getAltCount() != 0 ){
			//if the random allele picked is the reference, increase the reference counter
			if(arr->ref==sampleRandomRefAltAllele(arr->ref,
							      arr->alt,
							      arr->vectorAlleles->at(j).getRefCount(), 
							      arr->vectorAlleles->at(j).getAltCount())){
			    refCountg1+=arr->vectorAlleles->at(j).getRefCount();
			}else{//otherwise the alt
			    altCountg1+=arr->vectorAlleles->at(j).getAltCount();
			}
		    }else{
			refCountg1+=arr->vectorAlleles->at(j).getRefCount();
			altCountg1+=arr->vectorAlleles->at(j).getAltCount();
		    }

		}
	    }
		    
	    for(unsigned k=0;k<g2i.size();k++){
		if(j==g2i[k]){
		    //skip undefined sites 
		    if(arr->vectorAlleles->at(j).getRefCount() == 0 && 
		       arr->vectorAlleles->at(j).getAltCount() == 0 )
			goto nextiterationshare;			       


		    if(arr->vectorAlleles->at(j).getRefCount() != 0 && //hetero, increase one at random
		       arr->vectorAlleles->at(j).getAltCount() != 0 ){
			//if the random allele picked is the reference, increase the reference counter
			if(arr->ref==sampleRandomRefAltAllele(arr->ref,
							      arr->alt,
							      arr->vectorAlleles->at(j).getRefCount(), 
							      arr->vectorAlleles->at(j).getAltCount())){
			    refCountg2+=arr->vectorAlleles->at(j).getRefCount();
			}else{//otherwise the alt
			    altCountg2+=arr->vectorAlleles->at(j).getAltCount();
			}
		    }else{
			refCountg2+=arr->vectorAlleles->at(j).getRefCount();
			altCountg2+=arr->vectorAlleles->at(j).getAltCount();
		    }
		}
	    }
		    

	}
		
	//print data row

	if (refCountg1 != 0 ){ //if first group has ref
	    if(refCountg2 == 0)  // we must observe ref in the second
		goto nextiterationshare;
	}

	if (refCountg2 != 0 ){ //if second group has ref
	    if(refCountg1 == 0)  // we must observe ref in the first
		goto nextiterationshare;
	}

	if (altCountg1 != 0 ){ //if first group has alt
	    if(altCountg2 == 0) // we must observe alt in the second
		goto nextiterationshare;
	}

	if (altCountg2 != 0 ){ //if second group has alt
	    if(altCountg1 == 0) // we must observe alt in the first
		goto nextiterationshare;
	}



	if(!gw->writeAlleleRecord(arr)){
	    cerr<<"GlacSharing: error record "<<*arr<<endl;
	    exit(1);
	}
	keptRecords++;
		
    nextiterationshare:	    
	continue;

    }


    delete(gw);


    cerr<<"Program sharing wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

