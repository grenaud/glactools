
#include "ACF2GPHOCS.h"


using namespace std;


ACF2GPHOCS::ACF2GPHOCS(){

}

ACF2GPHOCS::~ACF2GPHOCS(){

}

string ACF2GPHOCS::usage() const{
    string usage=string("")+"glactools acf2ghocs  [options]  <ACF file> <bedfile>\n\n"
			"\tThis program produces gPhocs sequence input based on the ACF file and from the regions defined in the bedfile.\n"+
                	"\tOptions:\n"+
                  	"\t\t--allowCpg\t\tOnly allow transversions (Default "+boolStringify(allowCpg)+" )\n"+
           
                	"\tWARNING: Use preferably the wrapper acf2gphocsWrapper.pl because it computes\n"+
			"\tthe bed regions the puts the # of records in the header\n"+
	"";
    return usage;
}

int ACF2GPHOCS::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {

    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    for(int i=1;i<(argc-2);i++){ 
        
        if( string(argv[i]) == "--allowCpg"  ){
            allowCpg=true;
            continue;
        }

        cerr<<"Error unknown option "<<argv[i]<<endl;
        exit(1);
    }




    string bedFileRegions = string(argv[argc-1]);
    map< string, vector<GenomicRange> * > * bedRegionsToFilter;
    bedRegionsToFilter = readBEDSortedfile(bedFileRegions);

    map< string, unsigned int > * coordinateOfVec=new map< string, unsigned int >();
    for(map<string,vector<GenomicRange> * >::iterator it = bedRegionsToFilter->begin(); 
	it != bedRegionsToFilter->end(); 
	++it) {
	//cout<<it->first<<endl;
	coordinateOfVec->insert( pair<string ,unsigned int>(it->first,0) );
    }

    GlacParser gp  (argv[argc-2]);
    string chrName="-1";//current chr name 


    unsigned int previousCoordinate=0;
    vector<GenomicRange> * currentGr=0;
    unsigned int currentIndex=0;
    AlleleRecords * dataRow;
    bool chrFoundInBed=false;



    string         chrPrinted                = "-1";
    unsigned int   previousCoordinatePrinted =    0;
    bool           inBlockOfSequence         =false;
    unsigned int   locusNumber               =    0;
    vector<string> sequencesToPrint (gp.getPopulationsNames()->size()-1,"");//do not produce anc
    unsigned int   keptRecords               =    0;
    unsigned int   totalRecords              =    0;

    while(gp.hasData()){
	dataRow = gp.getData();


	if(dataRow->chr != chrName){//new chr
	    if(inBlockOfSequence){
		//TODO flush block

		inBlockOfSequence=false;
	    }

	    if(bedRegionsToFilter->find(dataRow->chr) == bedRegionsToFilter->end() ){
		chrFoundInBed=false;
		currentGr    = 0;		       
	    }else{
		chrFoundInBed=true;
		currentGr    = bedRegionsToFilter->at(dataRow->chr) ;
		currentIndex = coordinateOfVec->at(   dataRow->chr) ;
		if(currentIndex!=0){
		    cerr<<"There seems to be a mix of chromosomes in the mistar file, needs to be sorted chr: "<<chrName<<endl;
		    return 1;
		}
	    }
	    previousCoordinate = dataRow->coordinate;
	    chrName            = dataRow->chr;

	}else{// not new chr

	    if(previousCoordinate == dataRow->coordinate){
		cerr<<"WARNING: There seems to be a unsorted coordinate in the mistar file, needs to be sorted coordinate: "<<chrName<<":"<<previousCoordinate<<endl;
		//return 1;
		continue;
	    }

	    if(previousCoordinate >  dataRow->coordinate){
		cerr<<"ERROR: There seems to be a unsorted coordinate in the mistar file, needs to be sorted coordinate: "<<chrName<<":"<<previousCoordinate<<endl;
		return 1;
	    }
	}

	if(!chrFoundInBed )
	    goto nextmistarrecord;

	

	while(1){
	    if(currentIndex == currentGr->size())//end of bed file
		goto nextmistarrecord;
	    //ignore
	    //        |---------|
	    // *     
	    if( dataRow->coordinate < currentGr->at(currentIndex).getStartCoord() ){//ignore
		// cout<<"case 1"<<endl;
		if(inBlockOfSequence){
		    //flushing
		    cout<<"locus"<<++locusNumber<<" "<<(gp.getPopulationsNames()->size()-1)<<" "<<sequencesToPrint[0].size()<<endl;//minus anc
		    
		    for(unsigned int p=0;p<gp.getPopulationsNames()->size();p++){
			if(p==1) continue;//no anc
			if(p==0) 
			    cout<<gp.getPopulationsNames()->at(p)<<"\t"<<sequencesToPrint[p  ]<<endl;
			else
			    cout<<gp.getPopulationsNames()->at(p)<<"\t"<<sequencesToPrint[p-1]<<endl;
		    }
		    cout<<endl;
		    inBlockOfSequence=false;
		}


		goto nextmistarrecord;
	    }

	    //print
	    //        |---------|
	    //            *     
	    if(dataRow->coordinate >= currentGr->at(currentIndex).getStartCoord() &&
	       dataRow->coordinate <= currentGr->at(currentIndex).getEndCoord() ){
#ifdef DEBUG
		cerr<<"case 2"<<endl;
		cerr<<inBlockOfSequence<<" "<<(*dataRow)<<endl;
		cerr<<dataRow->chr<<"\t"<<dataRow->coordinate<<"\t"<<currentIndex<<endl;
#endif
		char refB=dataRow->ref;
		char altB=dataRow->alt;
		char hetB='N';//UIPAC base for het sites

		if(altB!='N')
		    hetB=dinucleotide2uipac(refB,altB);


		
		if(!inBlockOfSequence){//first time in block
		    //init
		    sequencesToPrint = vector<string>  (gp.getPopulationsNames()->size()-1,"");//init, do not produce anc and root
		    inBlockOfSequence=true;

		   
		    chrPrinted                = dataRow->chr        ;
		    previousCoordinatePrinted = dataRow->coordinate ;
		    
		}else{
		    //we were already in a block, inblockofsequence is true

		    if( (dataRow->chr        != chrPrinted) ){
			//end of a block but overshot to next chr

			//flush current locus
			cout<<"locus"<<++locusNumber<<" "<<(gp.getPopulationsNames()->size()-1)<<" "<<sequencesToPrint[0].size()<<endl;


			for(unsigned int p=0;p<gp.getPopulationsNames()->size();p++){
			    if(p==1) continue;//no anc
			    if(p==0) 
				cout<<gp.getPopulationsNames()->at(p)<<"\t"<<sequencesToPrint[p  ]<<endl;
			    else
				cout<<gp.getPopulationsNames()->at(p)<<"\t"<<sequencesToPrint[p-1]<<endl;
			    
			}
			cout<<endl;

			//re-init
			sequencesToPrint = vector<string>  (gp.getPopulationsNames()->size()-1,"");//init, do not produce anc and root

			//leave inBlockOfSequence flag as is
		    }

		    //check if we need print additional records
		    if( (dataRow->chr        == chrPrinted) &&
		 	(dataRow->coordinate != (previousCoordinatePrinted+1) ) ){

			//cout<<"MISSING"<<chrPrinted<<":"<<previousCoordinatePrinted<<endl;

			if( dataRow->coordinate <= (previousCoordinatePrinted+1)){
			    cerr<<"ERROR: mistar file does not appear to be sorted, died at coordinate "<<dataRow->chr<<":"<<dataRow->coordinate<<endl;
			    exit(1);
			}
			
			//print Ns in between 
			for(unsigned int i=(previousCoordinatePrinted+1);i<=(dataRow->coordinate-1);i++){
			    //cout<<"ADD MISSING"<<chrPrinted<<":"<<i<<endl;
			    for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
				if(j==1) continue;
				if(j==0)
				    sequencesToPrint[j  ] += "N";
				else
				    sequencesToPrint[j-1] += "N";				
			    }
			}
			//end of printing intervening Ns

		    }
		    
		    			    		       			
		} //end else already in block


		//add current bases
#ifdef DEBUG
		cerr<<"add "<<endl;
		cerr<<inBlockOfSequence<<" "<<(*dataRow)<<endl;
#endif
		previousCoordinatePrinted = dataRow->coordinate;
		chrPrinted                = dataRow->chr;

		bool isCpg=false;
		if(!allowCpg){
		    for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
			isCpg = isCpg || dataRow->vectorAlleles->at(j).getIsCpg();
		    }
		}

		if(isCpg){//mask CpG
		    
		    for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
			if(j==1) continue;		       
			if(j==0){//exception for the root, can be 1,0 or 0,1			    
			    sequencesToPrint[j ] += "N";			    
			}else{
			    sequencesToPrint[j-1] += "N";
			}
			continue;			
		    }
		    
		}else{
		
		    for(unsigned j=0;j<dataRow->vectorAlleles->size();j++){
			if(j==1) continue;

			//undefined site
			if( (dataRow->vectorAlleles->at(j).getRefCount() + dataRow->vectorAlleles->at(j).getAltCount()) > 2) {
			    cerr<<"ERROR: population "<<j<<" has more than 2 alleles at coordinate "<<dataRow->chr<<":"<<dataRow->coordinate<<endl;
			    exit(1);
			}

			if( (dataRow->vectorAlleles->at(j).getRefCount() + dataRow->vectorAlleles->at(j).getAltCount()) != 2) {
			    if(j==0){//exception for the root, can be 1,0 or 0,1
				if( (dataRow->vectorAlleles->at(j).getRefCount() == 1) && (dataRow->vectorAlleles->at(j).getAltCount() == 0) ){ 
				    sequencesToPrint[j ] += refB;
				    continue;
				}
			    
				if( (dataRow->vectorAlleles->at(j).getRefCount() == 0) && (dataRow->vectorAlleles->at(j).getAltCount() == 1) ){
				    sequencesToPrint[j ] += altB;
				    continue;
				}
			    
				sequencesToPrint[j ] += "N";			    
			    }else{
				sequencesToPrint[j-1] += "N";
			    }
			    continue;
			}
			//cout<<"add2 "<<sequencesToPrint[j-2]<<" "<<refB<<endl;
		    
			//homo ref
			if( (dataRow->vectorAlleles->at(j).getRefCount() == 2) && 
			    (dataRow->vectorAlleles->at(j).getAltCount() == 0) ){
			    if(j==0)
				sequencesToPrint[j  ] += refB;
			    else
				sequencesToPrint[j-1] += refB;
			    continue;
			}
			//cout<<"add3 "<<sequencesToPrint[j-2]<<" "<<hetB<<endl;

			//het ref+alt
			if( (dataRow->vectorAlleles->at(j).getRefCount() == 1) && 
			    (dataRow->vectorAlleles->at(j).getAltCount() == 1) ){
			    if(j==0)
				sequencesToPrint[j  ] += hetB;
			    else
				sequencesToPrint[j-1] += hetB;
			    continue;
			}
			//cout<<"add4 "<<sequencesToPrint[j-2]<<" "<<altB<<endl;
			//homo alt
			if( (dataRow->vectorAlleles->at(j).getRefCount() == 0) && 
			    (dataRow->vectorAlleles->at(j).getAltCount() == 2) ){
			    if(j==0)
				sequencesToPrint[j  ] += altB;
			    else
				sequencesToPrint[j-1] += altB;
			    continue;
			}

		    }//end each data row
		
		}
		    

	 

		
		//
		keptRecords++;				
		goto nextmistarrecord;
	    }//end in bed record

	    //we are running behind in the bed file
	    //        |---------|
	    //                      *     
	    if(  dataRow->coordinate > currentGr->at(currentIndex).getEndCoord()){
		//cout<<"case 3\t"<<currentIndex<<"\t"<<currentGr->size()<<endl;
		
		if(inBlockOfSequence){

		    cout<<"locus"<<++locusNumber<<" "<<(gp.getPopulationsNames()->size()-1)<<" "<<sequencesToPrint[0].size()<<endl;
		    for(unsigned int p=0;p<gp.getPopulationsNames()->size();p++){
			if(p==1) continue;
			
			if(p==0)
			    cout<<gp.getPopulationsNames()->at(p)<<"\t"<<sequencesToPrint[p  ]<<endl;
			else
			    cout<<gp.getPopulationsNames()->at(p)<<"\t"<<sequencesToPrint[p-1]<<endl;
		    }
		    cout<<endl;

		    inBlockOfSequence=false;
		}

		if(currentIndex<currentGr->size()){//we move to next iteration			    
		    currentIndex++;
		}else{
		    goto nextmistarrecord; //we have reached the end of the vector, do nothing until next chr
		}
	    }
	}
    nextmistarrecord:
	// cout<<dataRow->coordinate<<endl;		    

	totalRecords++;
    }

    cerr<<"Program "<<argv[0]<<" terminated gracefully, wrote "<<keptRecords<<" records out of  "<<totalRecords<<" records"<<endl;

    return 0;
}

