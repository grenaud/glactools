
#include "Glac2BED.h"


using namespace std;


Glac2BED::Glac2BED(){

}

Glac2BED::~Glac2BED(){

}

string Glac2BED::usage() const{
    string usage=string("glactools")+" glac2bed  <ACF or GLF file>"+
	"\nThis program will print the records regions covered in the ACF/GLF file as bed file\n\n";
    return usage;
}

int Glac2BED::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    string glacfile  = string(argv[argc-1]);
    
    GlacParser gp (glacfile);
    
	



    string chrName="-1";
    unsigned int previousCoordinate=0;
    unsigned int startBed=0;
    unsigned int endBed=0;
    

    AlleleRecords * dataRow;

    while(gp.hasData()){
	
	dataRow = gp.getData();

	    
	
	if(dataRow->chr != chrName){//new chr

	    if(chrName != "-1"){
		cout<<chrName<<"\t"<<startBed<<"\t"<<endBed<<endl;
	    }

	    previousCoordinate = dataRow->coordinate;
	    chrName            = dataRow->chr;
	    startBed = dataRow->coordinate-1;
	    endBed   = dataRow->coordinate;


	}else{
	 
	    if(previousCoordinate == dataRow->coordinate){
		cerr<<"WARNING: There seems to be redundant records in the ACF file, see coordinate: "<<chrName<<":"<<previousCoordinate<<endl;
		//return 1;
	    }

	    if(previousCoordinate > dataRow->coordinate){
		cerr<<"ERROR: There seems to be a unsorted coordinate in the ACF file, needs to be sorted coordinate: "<<chrName<<":"<<previousCoordinate<<endl;
		return 1;
	    }

	    if( (previousCoordinate+1) == dataRow->coordinate){
		endBed++;
	    }else{
		cout<<chrName<<"\t"<<startBed<<"\t"<<endBed<<endl;
		startBed= dataRow->coordinate-1;
		endBed  = dataRow->coordinate;
	    }

	    previousCoordinate = dataRow->coordinate;

	}
    }	
    cout<<chrName<<"\t"<<startBed<<"\t"<<endBed<<endl;

	
    return 0;
}

