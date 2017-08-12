
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
    for(int i=1;i<(argc-2);i++){ 
	
	if( string(argv[i]) == "--noanc"){
	    printRoot=false;
	    continue;
	}

	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }

    string glacfile  = string(argv[argc-2]);
    string bedFile   = string(argv[argc-1])+".bed";
    string bimFile   = string(argv[argc-1])+".bim";
    string famFile   = string(argv[argc-1])+".fam";

    GlacParser gp (glacfile);
    
    if(!gp.isACFormat()){
	cerr<<"ACF2EIGENSTRAT: Error the file "<<glacfile<<" should be in ACF format"<<endl;
        return 1;	
    }

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

	if(counter!=0 && (counter%100000)==0){
	    cerr<<"ACF2BPLINK: at "<<record->chr<<"\t"<<record->coordinate<<endl;
	}
	bimFileS<<record->chr<<"\t"<<"snp#"<<(counter++)<<"\t"<<stringify(double(record->coordinate)/double(1000000))<<"\t"<<stringify(record->coordinate)<<"\t"<<record->ref<<"\t"<<record->alt<<endl;
	
	unsigned int firstIndex=0;
	if(!printRoot)
	    firstIndex=2;
	
	
	char byteToWrite=0;
	int storedInByte=0;
	int dataToWrite=0;
	for(unsigned int i_=firstIndex;i_<record->vectorAlleles->size();i_++){	    
	    unsigned int i=i_-firstIndex;
	    char   toStore = record->vectorAlleles->at(i_).printPLINK();	   
	    // cout<<"i "<<i<<" "<<(i%4)<<" "<<record->vectorAlleles->at(i_)<<endl;
	    // cout<<"1: "<<var2binary(toStore)<<endl;
	    if( (i%4) == 3){//to write
		byteToWrite |= (  toStore<< ((i%4)*2) );		
		//cout<<"3: "<<var2binary(byteToWrite)<<endl;
		bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		
		byteToWrite  = 0;
		storedInByte = 0;
		dataToWrite=0;
	    }else{
		dataToWrite++;
		byteToWrite |= (  toStore<< ((i%4)*2) );
		//cout<<"1: "<<var2binary(byteToWrite)<<endl;
		storedInByte++;
	    }
	} 

	if(dataToWrite != 0){
	    int k=3;
	    for(int i=0;i<(4-dataToWrite);i++){
		char zeroC= 0; //NULL
		byteToWrite |= (  zeroC << ((k%4)*2) );
		k--;
	    }
	    //cout<<"EX: "<<var2binary(byteToWrite)<<endl;  
	    bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		
	}
	//     }
	//     // cout<<"1: "<<var2binary(byteToWrite)<<endl;  
	//     bedFileS.write( (char *)&byteToWrite, sizeof(byteToWrite));		

	// }

	// if(storedInByte!=0){
	// }
	
	//bedFileS<<endl;

    }
    bedFileS.close();
    bimFileS.close();
    
    cerr<<"Program "<<argv[0]<<" terminated gracefully, wrote "<<counter<<" sites"<<endl;
	
    return 0;
}

