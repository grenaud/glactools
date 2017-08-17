
#include "GlacReplaceAncestor.h"


using namespace std;


GlacReplaceAncestor::GlacReplaceAncestor(){

}

GlacReplaceAncestor::~GlacReplaceAncestor(){

}

string GlacReplaceAncestor::usage() const{
    string usage=string("glactools")+" replaceanc [options] <GLAC file1>  <GLAC file2>"+
	"\nThis program will print the first GLAC file but with the ancestral information from the second one to STDOUT.\n"+
	"\n"+
	"Options:\n"+
	"\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n";
    return usage;
}

int GlacReplaceAncestor::run(int argc, char *argv[]){

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
	       

    if(lastOpt != (argc-2)){
	cerr<<"GlacReplaceAncestor: The last arguments are the <ACF file> "<<endl;
	return 1;		
    }


    string g1 = string(argv[lastOpt+0]);
    string g2 = string(argv[lastOpt+1]);



    GlacParser gp1 (g1);

    GlacParser gp2 (g2);
    // cerr<<"gp1 "<<gp1.getSizePops()<<endl;
    // cerr<<"gp2 "<<gp2.getSizePops()<<endl;

    AlleleRecords * arr;
    if(gp1.getHeaderSQ("") != gp2.getHeaderSQ("") ){
	cerr<<"GlacReplaceAncestor: The SQ fields differ in the header, are they from the same reference?"<<endl;
	return 1;
    }

    if(gp1.isGLFormat()){
	if(!gp2.isGLFormat()){	
	    cerr<<"GlacReplaceAncestor: ERROR: The first file is GLF but the second is ACF, "<<endl;
	    return 1;
	}
    }
    

    if(!gp1.isGLFormat()){
	if(gp2.isGLFormat()){	
	    cerr<<"GlacReplaceAncestor: ERROR: The first file is GLF but the second is ACF, "<<endl;
	    return 1;
	}
    }


    GlacWriter * gw = new GlacWriter(gp1.getSizePops(),
				     gp1.isGLFormat(),
				     gp1.isGLFormat()?1:2,
				     1,//compression threads
				     uncompressed);

    stringstream newheader;
    if(gp1.isGLFormat())
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
    newheader<<"#REPLACEANC:"<<"\n";

    newheader<<"#REPLACEANC#"<<1<<"\n";
    newheader<<gp1.getHeaderNoSQNoDefline("#\t")<<"\n";
    newheader<<"#REPLACEANC#"<<2<<"\n";
    newheader<<gp2.getHeaderNoSQNoDefline("#\t")<<"\n";

    newheader<<gp1.getHeaderSQ("")<<"\n";
    newheader<<gp1.getDefline()<<"\n";

    
    if(!gw->writeHeader(newheader.str())){
	cerr<<"GlacReplaceAncestor: error writing header "<<endl;
	return 1;
    }
	





    //cout<<vectorToString( *(mp1.getPopulationsNames())," ")<<" "<<vectorToString( *(gp2.getPopulationsNames())," ")<<endl;
    bool hasData1 = gp1.hasData();
    bool hasData2 = gp2.hasData();

    unsigned int nonPop1 = gp1.getPopulationsNames()->size();
    //unsigned int nonPop2 = gp2.getPopulationsNames()->size();
    

    //printing first
    // for(unsigned int i=2;i<nonPop1;i++){
    // 	cout<<gp1.getPopulationsNames()->at(i);
    // 	if( i!= (nonPop1-1) ){
    // 	    cout<<"\t";
    // 	}
    // }
    // cout<<endl;;
    // for(unsigned int i=1;i<nonPop2;i++){
    // 	cout<<gp2.getPopulationsNames()->at(i);
    // 	if(i!=(nonPop2-1)){
    // 	    cout<<"\t";
    // 	}
    // }


    if(!hasData1){
	cerr<<"File #1 does not have any data, exiting"<<endl;
	return 1;
    }
    if(!hasData2){
	cerr<<"File #2 does not have any data, exiting"<<endl;
	return 1;
    }

    AlleleRecords * record1 = gp1.getData();
    AlleleRecords * record2 = gp2.getData();

    if(record1->chri != record2->chri ){
	cerr<<"Warning: Initial chromosomes differ between "<<*record1<<" and "<<*record2<<", make sure it is the correct species/genome build"<<endl;
	//return 1;
    }

    // unsigned int coordCurrent=min( record1->coordinate,
    // 				   record2->coordinate );
    bool stayLoop=true;


    while(stayLoop){
	//cerr<<hasData1<<" "<<hasData2<<" "<<	   record1->chri <<" "<<record2->chri  <<" "<<record1->coordinate <<" "<<record2->coordinate<<endl;

	if(!hasData1  ){
	    stayLoop=false;
	    break;
	}


	//	return 1;

	//second one has data
	if(hasData2 &&
	   record1->chri       == record2->chri        &&
	   record1->coordinate == record2->coordinate ){
	    //cerr<<"SYNC "<<endl;
	    // return 1;
	    if(record1->chri != record2->chri ){
		cerr<<"Chromosomes differ between "<<(*record1)<<" and "<<(*record2)<<endl;
		return 1;
	    }

	    if(record1->ref != record2->ref ){
		cerr<<"The reference allele differs between "<<(*record1)<<" and "<<(*record2)<<endl;
		return 1;
	    }


	    
	    char newAlt='N';
	    // int rootLocated=0; //o=unknown,1=record1,2=record2,3=both should agree
	    // return 1;

	    if(record1->alt == record2->alt ){ //agree
		// rootLocated=3;
		newAlt=record1->alt;
		goto printnewrecord;
	    }
	    

	    if( record1->alt == 'N' && isResolvedDNA(record2->alt) ){
		// rootLocated=2;
		newAlt=record2->alt;
		goto printnewrecord;	
	    }

	    if(isResolvedDNA(record1->alt) && record2->alt == 'N' ){
		// rootLocated=1;
		newAlt=record1->alt;
		goto printnewrecord;	
	    }


	    //check for diff alt allele
	    if(record1->alt != record2->alt &&
	       isResolvedDNA(record1->alt)  &&
	       isResolvedDNA(record2->alt) ){ 

		//check if the alternative is due to the root or other populations
		//int sumAlt=0;
		bool oneHasAlt=false;
		for(unsigned int i=0;i<2;i++){	
	    
		    if(gp2.isGLFormat()){
			oneHasAlt=oneHasAlt || (*record2->vectorGLs)[i].hasAlt();
		    }else{
			oneHasAlt=oneHasAlt || (*record2->vectorAlleles)[i].hasAlt();
		    }
		    //		    SingleGL::hasAlt
		    //sumAlt+=(*record2->vectorAlleles)[i].getAltCount();
		}
		
		//if(sumAlt!=0) //the alt is either the root or anc
		if(oneHasAlt)
		    goto seekdata;//cannot reconcile alternative alleles

		//else safe ignore the second alternative since the second does not have the alt and pick the first as the alt
		newAlt=record1->alt;
		goto printnewrecord;//cannot reconcile alternative alleles
				    
		
	    }else{
		cerr<<"GlacReplaceAncestor: WARNING: Wrong state between "<<(*record1)<<" and "<<(*record2)<<endl;
		//return 1;
	    }

	    
	printnewrecord:
	    //cerr<<"printnewrecord"<<endl;
	    AlleleRecords arw;
	    arw.chri        = record1->chri;
	    arw.coordinate  = record1->coordinate;
	    arw.ref         = record1->ref;
	    arw.alt         = newAlt;

	    if(gp1.isGLFormat()){
		arw.vectorAlleles = 0;
		arw.vectorGLs     = new vector<SingleGL>     (  );
		for(unsigned int i=0;i<2;i++)
		    arw.vectorGLs->push_back( record2->vectorGLs->at(i) );		
		for(unsigned int i=2;i<record1->vectorGLs->size();i++)
		    arw.vectorGLs->push_back( record1->vectorGLs->at(i) );		
	    }else{
		arw.vectorAlleles = new vector<SingleAllele> (  );
		arw.vectorGLs     = 0;
		for(unsigned int i=0;i<2;i++)
		    arw.vectorAlleles->push_back( record2->vectorAlleles->at(i) );		    
		for(unsigned int i=2;i<record1->vectorAlleles->size();i++)
		    arw.vectorAlleles->push_back( record1->vectorAlleles->at(i) );		
	    }
	    
	    if(!gw->writeAlleleRecord(&arw)){
		cerr<<"GlacReplaceAncestor: error writing record "<<arw<<endl;
		return 1;
	    }
	    	    
	}
	//closes
	// hasData2 &&
	//    record1->chr        == record2->chr        &&
	//    record1->coordinate == record2->coordinate )
	else{//no record in the second one or not synched
	    
	    //check second record
 	    if(hasData2){
		//file 2 is behind, need to increase the 
		int chrcmp = compare2ChrsU(record1->chri,record2->chri);
		//chr1, vecAlleleRecords[i]->chr);

		//if( (record1->chr         >  record2->chr) ||
		if(  (chrcmp == 1) ||
		    ((record1->chri        == record2->chri) &&  ( record1->coordinate > record2->coordinate ) )
		){
		    hasData2 = gp2.hasData();
		    if(hasData2){
			record2 = gp2.getData();
		    }
		    //BEGIN TODELETE
		    // if( (record2->coordinate%100000) == 0){
		    // 	cerr<<"file 2 is behind "<<record2->chr<<":"<<record2->coordinate<<endl;			
		    // }
		    //END TODELETE
		    continue;//next iteration
		}
		    
 	    }

	    AlleleRecords arw;
	    arw.chri        = record1->chri;
	    arw.coordinate  = record1->coordinate;
	    arw.ref         = record1->ref;
	    //arw.alt         = newAlt;

	    // cout<<record1->chr<<"\t"<<record1->coordinate<<"\t"; //<<record2->coordinate<<"\t";
	    // cout<<record1->ref<<",";
	    // // cerr<<record1->chr<<"\t"<<record1->coordinate<<"\t"; //<<record2->coordinate<<"\t";
	    // // cerr<<record1->ref<<",";
	    
	    // return 1;
	    //int sumAlt=0;
	    bool oneHasAlt=false;
	    for(unsigned int i=2;i<record1->vectorAlleles->size();i++){
		//sumAlt+=(*record1->vectorAlleles)[i].getAltCount();
		if(gp1.isGLFormat()){
		    oneHasAlt=oneHasAlt || (*record1->vectorGLs)[i].hasAlt();
		}else{
		    oneHasAlt=oneHasAlt || (*record1->vectorAlleles)[i].hasAlt();
		}		
	    }
	    
	    //if(onsumAlt==0){ //the only alternative allele was due to the root
	    if(!oneHasAlt){ //the only alternative allele was due to the root
		//cout<<"N\t";		
		arw.alt = 'N';
	    }else{
		//cout<<record1->alt<<"\t";
		arw.alt = record1->alt;
	    }

	    if(gp1.isGLFormat()){
		arw.vectorAlleles = 0;
		arw.vectorGLs     = new vector<SingleGL>     (  );
		for(unsigned int i=0;i<2;i++){
		    SingleGL sgl;
		    arw.vectorGLs->push_back( sgl );
		}
	    }else{
		arw.vectorAlleles = new vector<SingleAllele> (  );
		arw.vectorGLs     = 0;
		for(unsigned int i=0;i<2;i++){
		    SingleAllele sal;
		    arw.vectorAlleles->push_back( sal );
		}

	    }
	    

	    if(gp1.isGLFormat()){
		for(unsigned int i=2;i<record1->vectorGLs->size();i++)
		    arw.vectorGLs->push_back( record1->vectorGLs->at(i) );		
	    }else{
		for(unsigned int i=2;i<record1->vectorAlleles->size();i++)
		    arw.vectorAlleles->push_back( record1->vectorAlleles->at(i) );		
	    }

	    if(!gw->writeAlleleRecord(&arw)){
		cerr<<"GlacReplaceAncestor: error writing record "<<arw<<endl;
		return 1;
	    }

	}

    seekdata:	
	if(hasData1){
	    hasData1 = gp1.hasData();
	    if(hasData1){
		record1 = gp1.getData();
	    }
	}



    }//end stayloop




    delete(gw);


    cerr<<"Program finished gracefully"<<endl;

    return 0;
}

