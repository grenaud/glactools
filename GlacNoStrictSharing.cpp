
#include "GlacNoStrictSharing.h"


using namespace std;


GlacNoStrictSharing::GlacNoStrictSharing(){

}

GlacNoStrictSharing::~GlacNoStrictSharing(){

}

string GlacNoStrictSharing::usage() const{
 string usage=string("glactools")+" snosharing  <ACF file>  <group 1> <group 2>"+

	"\nThis will filter sites where individuals in population group 1 do not share an allele with individual in population group 2.\n"+
	"The individuals/populations in <group 1> and <group 2> are comma-separated\n"
	"In other words:\n"+
	"If the individuals in group 1 have the reference allele, it must be must not be observed in group2\n"+
	"If the individuals in group 1 have the alternative allele, it must be not be observed in group2\n"+
	"Example: glactools snosharing input.acf.gz \"CEU,TSI\" \"CHB,JPT\"\n"+
	"It requires that the allele count for every individual for both groups be non-zero.\n"+
	"Unlike the \"sharing\", a random allele will not be picked for heterozygous (or variable in a population) positions so it will \n"+
	"effectively throw out every variable site in <group 1> <group 2>\n"+
	"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n";

    return usage;
}

int GlacNoStrictSharing::run(int argc, char *argv[]){

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
	cerr<<"GlacNoStrictSharing: The last arguments are the <ACF file> "<<endl;
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
	cerr<<"GlacNoStrictSharing: This function is only defined for  <ACF file> "<<endl;
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
    newheader<<"#STRICTNOSHARING:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacNoStrictSharing: error writing header "<<endl;
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
			goto nextiterationstrictnoshare;			       
			    
		    if(arr->vectorAlleles->at(j).getRefCount() != 0 && //hetero/variable skip
		       arr->vectorAlleles->at(j).getAltCount() != 0 ){
			goto nextiterationstrictnoshare;			       

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
			goto nextiterationstrictnoshare;			       

		    if(arr->vectorAlleles->at(j).getRefCount() != 0 && //hetero/variable skip
		       arr->vectorAlleles->at(j).getAltCount() != 0 ){
			goto nextiterationstrictnoshare;			       
		    }else{
			refCountg2+=arr->vectorAlleles->at(j).getRefCount();
			altCountg2+=arr->vectorAlleles->at(j).getAltCount();
		    }


		}
	    }
		    

	}
		
	//print data row

	if (refCountg1 != 0 && refCountg2 != 0) // if they share reference
	    goto nextiterationstrictnoshare;

	if (altCountg1 != 0 && altCountg2 != 0) // if they share alternative
	    goto nextiterationstrictnoshare;


	if(!gw->writeAlleleRecord(arr)){
	    cerr<<"GlacSegsite: error record "<<*arr<<endl;
	    exit(1);
	}
	keptRecords++;
		
    nextiterationstrictnoshare:	    
	continue;

    }


    delete(gw);


    cerr<<"Program segsite wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

