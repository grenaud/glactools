
#include "ACF2BPLINK.h"


using namespace std;


ACF2BPLINK::ACF2BPLINK(){

}

ACF2BPLINK::~ACF2BPLINK(){

}

string ACF2BPLINK::usage() const{
    string usage=string("glactools")+" acf2bplink  [options] <ACF file> [out prefix]"+
	"\nThis program takes a ACF file and prints the genotype and SNP file as binary PLINK files\n\n"+
	"\tOptions\n"+			
       "\t\t"+"--noanc"+"\t"+"Do not print the root/anc (Default: "+boolStringify(printRoot)+" )\n\n\n";
    return usage;
}

int ACF2BPLINK::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc < 2  ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    //all but last 2
    for(int i=0;i<(argc-2);i++){ 

	if( string(argv[i]) == "--noanc"){
	    printRoot=false;
	}
	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    string glacfile  = string(argv[argc-2]);
    string bedFile   = string(argv[argc-1])+".bed";
    string bimFile   = string(argv[argc-1])+".bim";
    string famFile   = string(argv[argc-1])+".fam";

    GlacParser gp (glacfile);
    


    ofstream bedFileS;
    ofstream bimFileS; 
    ofstream famFileS; 
	
    bedFileS.open(bedFile.c_str(),   ios::out| ios::binary);
    bimFileS.open(bimFile.c_str(),   ios::out);
    famFileS.open(famFile.c_str(),   ios::out);

    
    if (!bedFileS.good()){  cerr << "Unable to open file "<<bedFile<<endl; return 1;   }
    if (!bimFileS.good()){  cerr << "Unable to open file "<<bimFile<<endl; return 1;   }
    if (!famFileS.good()){  cerr << "Unable to open file "<<famFile<<endl; return 1;   }
    

    AlleleRecords * record;
 
    unsigned int firstIndex=0;

    if(!printRoot)
	firstIndex=2;

    for(unsigned int i=firstIndex;i<gp.getPopulationsNames()->size();i++){
    	famFileS<<gp.getPopulationsNames()->at(i)<<"\t"<<gp.getPopulationsNames()->at(i)<<"\t"<<gp.getPopulationsNames()->at(i)<<"\t"<<gp.getPopulationsNames()->at(i)<<"\t1\t-9"<<endl;	
    } 
    famFileS.close();

    //magic number
    char c = 108;
    bedFileS.write( (char *)&c, sizeof(c));
    c = 27;
    bedFileS.write( (char *)&c, sizeof(c));

    //snp major
    c = 1;
    bedFileS.write( (char *)&c, sizeof(c));



    unsigned int counter=0;
    while(gp.hasData()){
	record = gp.getData();
	
	if(!isResolvedDNA(record->alt))
	    continue;

	bimFileS<<record->chr<<"\t"<<"snp#"<<(counter++)<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	unsigned int firstIndex=0;
	if(!printRoot)
	    firstIndex=2;
	

	char byteToWrite=0;
	int storedInByte=0;
	for(unsigned int i=firstIndex;i<record->vectorAlleles->size();i++){	    
	    char   toStore = record->vectorAlleles->at(i).printPLINK();	   
	    if( (i%4) == 3){//to write
		byteToWrite |= (  toStore<< ((i%4)*2) );		
		// cout<<"1: "<<var2binary(byteToWrite)<<endl;
		bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		
		byteToWrite  = 0;
		storedInByte = 0;
	    }else{
		byteToWrite |= (  toStore<< ((i%4)*2) );
		storedInByte++;
	    }
	} 

	if(storedInByte!=0){
	    for(int k=3;k>=storedInByte;k--){
		char zeroC= 3; //missing data
		byteToWrite |= (  zeroC << ((k%4)*2) );		
	    }
	    // cout<<"1: "<<var2binary(byteToWrite)<<endl;  
	    bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		
	}
	
	//bedFileS<<endl;

    }
    bedFileS.close();
    bimFileS.close();
    
    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;
	
    return 0;
}

