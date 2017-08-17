
#include "GlacClosest.h"


using namespace std;


GlacClosest::GlacClosest(){

}

GlacClosest::~GlacClosest(){

}

string GlacClosest::usage() const{
    string usage=string("glactools")+" closest  <ACF or GLF file>"+
	"\nThis program will print to stdout the distance to the closest site for each record\n\n";
    return usage;
}

int GlacClosest::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    string glacfile  = string(argv[argc-1]);
    
    GlacParser gp (glacfile);               

    AlleleRecords * dataRow;
    unsigned int totalRecords=0;
    unsigned int lastCoordinateM2=0;
    unsigned int lastCoordinateM1=0;
    unsigned int lastCoordinateM =0;
    uint16_t chri=UINT16_MAX;


    bool newChr1=true;
    bool newChr2=true;


    while(gp.hasData()){
	
	dataRow = gp.getData();


	if( chri != dataRow->chri){
	    //cerr<<"Chromosomes differ for line :"<<(*dataRow)<<endl;
	    newChr1=true;
	    newChr2=true;
	    lastCoordinateM2=0;
	    lastCoordinateM1=0;
	    lastCoordinateM =0;
	    //return 1;
	}

	if(newChr1){
	    lastCoordinateM2=dataRow->coordinate;
	    chri=dataRow->chri;
	    totalRecords++;
	    newChr1=false;

	    continue;
	}

	if(newChr2){
	    lastCoordinateM1=dataRow->coordinate;
	    if(lastCoordinateM1<lastCoordinateM2){
		cerr<<"Coordinate are not sorted :"<<(*dataRow)<<endl;
		return 1;
	    }

	    chri=dataRow->chri;
	    totalRecords++;
	    cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;
	    newChr2=false;

	    continue;
	}


	lastCoordinateM=dataRow->coordinate;

	if(lastCoordinateM<lastCoordinateM1){
	    cerr<<"Coordinate are not sorted :"<<(*dataRow)<<endl;
	    return 1;
	}

	if( (lastCoordinateM1-lastCoordinateM2) < (lastCoordinateM-lastCoordinateM1) ){
	    cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;
	}else{
	    cout<<(lastCoordinateM-lastCoordinateM1)<<endl;
	}
	
	//for next iteration
	lastCoordinateM2=lastCoordinateM1;
	lastCoordinateM1=lastCoordinateM;

	totalRecords++;
	
    }	

    cout<<(lastCoordinateM1-lastCoordinateM2)<<endl;
    
    cerr<<"Program "<<argv[0]<<" looked at "<<totalRecords<<" terminated gracefully"<<endl;


	
    return 0;
}

