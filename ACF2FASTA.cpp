
#include "ACF2FASTA.h"


using namespace std;


ACF2FASTA::ACF2FASTA(){

}

ACF2FASTA::~ACF2FASTA(){

}

string ACF2FASTA::usage() const{
    string usage=string("glactools")+" acf2fasta  [options] <ACF file>"+
	"\nThis program takes an ACF file and prints a FASTA file using the allele information\nwith one record per population. Each site generates one base pair.\n\n"+
	"\tOptions\n"+			
	"\t\t"+"--noanc"+"\t"+"Do not print the root/anc (Default: "+boolStringify(printRoot)+" )\n"
	"\t\t"+"--het"+"\t"+"Produce two fasta files containing the alleles for het sites (Default: "+boolStringify(produceTwo)+" )\n\n\n";
    return usage;
}

int ACF2FASTA::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {

    //all but last arg
    int lastOpt=1;

    for(int i=1;i<(argc);i++){ 

	if((string(argv[i]) == "-")  ){
            lastOpt=i;
            break;          
        }

        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }

	if( string(argv[i]) == "--noanc"){
	    printRoot=false;
	    continue;
	}

	if( string(argv[i]) == "--het"){
	    produceTwo=true;
	    continue;
	}

	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    string glacfile  = string(argv[lastOpt]);
    
    GlacParser gp (glacfile);
    
           






    AlleleRecords * record;
 
    unsigned int firstIndex=0;
    if(!printRoot)
	firstIndex=2;

    vector<string> deflines;
    vector<string> sequences;

    if(produceTwo){

	for(unsigned int i=firstIndex;i<gp.getPopulationsNames()->size();i++){
	    //indFileS<<gp.getPopulationsNames()->at(i)<<"\tU\t"<<gp.getPopulationsNames()->at(i)<<endl;	
	    deflines.push_back(gp.getPopulationsNames()->at(i)+"-1");
	    sequences.push_back("");
	    deflines.push_back(gp.getPopulationsNames()->at(i)+"-2");
	    sequences.push_back("");
	} 

    }else{
	for(unsigned int i=firstIndex;i<gp.getPopulationsNames()->size();i++){
	    //indFileS<<gp.getPopulationsNames()->at(i)<<"\tU\t"<<gp.getPopulationsNames()->at(i)<<endl;	
	    deflines.push_back(gp.getPopulationsNames()->at(i));
	    sequences.push_back("");
	} 
	// indFileS.close();
    }
    
    //unsigned int counter=0;
    while(gp.hasData()){
	record = gp.getData();
	// //	cout<<"ok"<<endl;    
	// cout<<*record<<endl;
	// // if(!isResolvedDNA(record->alt))
	// //     continue;

	// snpFileS<<"snp#"<<(counter++)<<"\t"<<record->chr<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	unsigned int firstIndex=0;
	if(!printRoot)
	    firstIndex=2;
	unsigned int indexVec=0;
	if(produceTwo){

	    for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
		char c = record->vectorAlleles->at(i).generateRandomAllele(record->ref,record->alt);//otherwise, a file will have the major allele and the other one the minor.
		char c2;
		if(record->vectorAlleles->at(i).isHeterozygous()){
		    sequences[ int(indexVec*2.0)   ]  += c;

		    if(record->ref == c)
			c2=record->alt;
		    else
			c2=record->ref;
		    
		    sequences[ int(indexVec*2.0)+1 ]  += c2;

		}else{
		    sequences[ int(indexVec*2.0)   ]  += c;
		    sequences[ int(indexVec*2.0)+1 ]  += c;		
		}
		    
		indexVec++;
	    }
	    
	}else{
	    for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){
		sequences[ indexVec++ ]  += record->vectorAlleles->at(i).generateRandomAlleleBias(record->ref,record->alt);
	    }
	}	
	// genoFileS<<endl;

    }
    // genoFileS.close();
    // snpFileS.close();

    for(unsigned int i=0;i<sequences.size();i++){
	cout<<">"<<deflines[i]<<endl;
	cout<<sequences[i]<<endl;
    }
    
    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

	
    return 0;
}

