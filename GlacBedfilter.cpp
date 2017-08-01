
#include "GlacBedfilter.h"


using namespace std;


GlacBedfilter::GlacBedfilter(){

}

GlacBedfilter::~GlacBedfilter(){

}

string GlacBedfilter::usage() const{
    string usage=string("glactools")+" bedfilter  <ACF or GLF file> <sorted BED file>"+
	"\nThis will keep only the positions in the bed file\n"+
"\n"+
	"Options:\n"+
    "\t"+"-u" +    "\t\t\t"+"Produce binary but uncompressed output (default: "+booleanAsString(uncompressed)+")\n";

    return usage;
}

int GlacBedfilter::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }
    int lastOpt=1;
    unsigned int firstPopInd=0;
		      

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

	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;
    }


    if(lastOpt != (argc-2)){
	cerr<<"The last 2 arguments are the <ACF/GLF file> <BED file> "<<endl;
	return 1;		
    }


    string fileglf        = string(argv[lastOpt]);
    string bedFileRegions = string(argv[lastOpt+1]);
    map< string, vector<GenomicRange> * > * bedRegionsToFilter;

    bedRegionsToFilter = readBEDSortedfile(bedFileRegions);
    map< string, uint32_t > * coordinateOfVec=new map< string, uint32_t >();
    for(map<string,vector<GenomicRange> * >::iterator it = bedRegionsToFilter->begin(); 
	it != bedRegionsToFilter->end(); 
	++it) {
	//cout<<it->first<<endl;
	coordinateOfVec->insert( pair<string ,uint32_t>(it->first,0) );
    }
    uint16_t chrNameIdx=-1;
    uint32_t  previousCoordinate=0;
    vector<GenomicRange> * currentGr=0;
    unsigned int currentIndex=0;
    bool chrFoundInBed=false;

    GlacParser gp (fileglf);
    AlleleRecords * arr;
    

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
    newheader<<"#BEDFILTER:"<<"\n";
    newheader<<gp.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<gp.getHeaderSQ("")<<"\n";
    newheader<<gp.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacBedfilter: error writing header "<<endl;
	exit(1);
    }
	

    uint64_t totalRecords=0;
    uint64_t keptRecords=0;
    
    
    while(gp.hasData()){
	arr = gp.getData();
	totalRecords++;


	if(arr->chri != chrNameIdx){//new chr
	    if(bedRegionsToFilter->find(gp.getChromosomeName(arr->chri)) == bedRegionsToFilter->end() ){
		chrFoundInBed=false;
		currentGr    = 0;		       
	    }else{
		chrFoundInBed=true;
		currentGr    = bedRegionsToFilter->at(arr->chr) ;
		currentIndex =   coordinateOfVec->at(arr->chr) ;
		if(currentIndex!=0){
		    cerr<<"There seems to be a mix of chromosomes in the mistar file, needs to be sorted chr: "<<chrNameIdx<<endl;
		    return 1;
		}
	    }
	    previousCoordinate = arr->coordinate;
	    chrNameIdx         = arr->chri;
	}else{
	    if(previousCoordinate >= arr->coordinate){
		cerr<<"There seems to be a unsorted coordinate in the mistar file, needs to be sorted coordinate: "<<previousCoordinate<<endl;
		return 1;
	    }
	}
	
	if(!chrFoundInBed )
	    goto nextacfrecord;
	
	
	
	while(1){
	    if(currentIndex == currentGr->size())//end of bed file
		goto nextacfrecord;
	    //ignore
	    //        |---------|
	    // *     
	    if( arr->coordinate < currentGr->at(currentIndex).getStartCoord() ){//ignore
		// cout<<"case 1"<<endl;
		goto nextacfrecord;
	    }
	    
	    //print
	    //        |---------|
	    //            *     
	    if(arr->coordinate >= currentGr->at(currentIndex).getStartCoord() &&
	       arr->coordinate <= currentGr->at(currentIndex).getEndCoord() ){
		// cout<<"case 2"<<endl;
		//cout<<(*arr)<<endl;
		if(!gw->writeAlleleRecord(arr)){
		    cerr<<"GlacBedfilter: error record "<<*arr<<endl;
		    exit(1);
		}
		keptRecords++;				
		goto nextacfrecord;
	    }
	    
	    //we are running behind in the bed file
	    //        |---------|
	    //                      *     
	    if(  arr->coordinate > currentGr->at(currentIndex).getEndCoord()){
		// cout<<"case 3\t"<<currentIndex<<"\t"<<currentGr->size()<<endl;
		
		if(currentIndex<currentGr->size()){//we move to next iteration			    
		    currentIndex++;
		}else{
		    goto nextacfrecord; //we have reached the end of the vector, do nothing until next chr
		}
	    }
	}
    nextacfrecord:
	// cout<<arr->coordinate<<endl;		    
	
	totalRecords++;
		// if(gp.isGLFormat()){
	//     for(unsigned j=firstPopInd;j<arr->vectorGLs->size();j++){
	// 	if( arr->vectorGLs->at(j).alleleCountNull() )
	// 	    goto nextiterationnoundef;		
	//     }
	// }else{
	//     for(unsigned j=firstPopInd;j<arr->vectorAlleles->size();j++){
	// 	if( arr->vectorAlleles->at(j).alleleCountNull() )
	// 	    goto nextiterationnoundef;
	//     }

	// }

	// if(!gw->writeAlleleRecord(arr)){
	//     cerr<<"GlacViewer: error record "<<*arr<<endl;
	//     exit(1);
	// }
	keptRecords++;

    nextiterationnoundef:	    
	continue;
	

    }


    delete(gw);


    cerr<<"Program undef wrote "<<keptRecords<<" out of "<<totalRecords<<" terminated gracefully"<<endl;

    return 0;
}

