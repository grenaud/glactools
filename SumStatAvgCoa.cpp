/*
 * SumStatAvgCoa
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SumStatAvgCoa.h"


SumStatAvgCoa::SumStatAvgCoa(){
    divergenceResults=NULL;
    populationNames =NULL;
}

SumStatAvgCoa::SumStatAvgCoa(const SumStatAvgCoa & other){
    // cout<<"copy"<<endl;
    numberOfPopulations = other.numberOfPopulations;
    // cout<<"copy "<<numberOfPopulations<<endl;

    populationNames = new vector<string>();

    for(unsigned int i=0;i<other.populationNames->size();i++){ 
	populationNames->push_back(other.populationNames->at(i));
	// cout<<"copyNames "<< other.populationNames->at(i) <<endl;
    }

    divergenceResults = new AvgCoaResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
	    divergenceResults[i][j]=other.divergenceResults[i][j];	    
	    // cout<<"copy "<<divergenceResults[i][j]<<endl;
	}
    }
    // cout<<"copy"<<endl;
}


SumStatAvgCoa::SumStatAvgCoa(const vector<string> * popNames){

    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    divergenceResults = new AvgCoaResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
    }
    populationNames->push_back("ref");
    
}

// SumStatAvgCoa::SumStatAvgCoa(const SumStatAvgCoa & other){

// }

SumStatAvgCoa::~SumStatAvgCoa(){
    // cout<<"destructor#1"<<endl;
    if(populationNames!=NULL)
	delete(populationNames);
    //cout<<"destructor#2"<<endl;
    if(divergenceResults!=NULL){
	for(unsigned int i=0;i<numberOfPopulations;i++){
	    //cout<<"destructor #"<<i<<endl;
	    delete [] divergenceResults[i];
	}
	delete [] divergenceResults;
    }
}

AvgCoaResult const * const *  SumStatAvgCoa::getAvgCoaResult() const{
    return divergenceResults;
}

string SumStatAvgCoa::printWithBootstraps(const   vector<SumStatAvgCoa *> * jackVec,  const string & dnaDistMode) const{ //all boostraps or jacknife
    // string toReturn;
    
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){

	//for each population, except the root/ancestral at index 0,1
	for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   

	    vector< const AvgCoaResult * > jacknivesToSend;
	    for(unsigned k=0;k<jackVec->size();k++){	       
		const AvgCoaResult * const* temsp= (jackVec->at(k)->getAvgCoaResult());
		//constAvgCoaResult
		//cout<<temsp[i][j]<<endl;
		//const AvgCoaResult * test1  = &temsp[i][j];
		//jacknivesToSend.push_back( &(  (jackVec->at(k)->getAvgCoaResult())[i][j]) );
		jacknivesToSend.push_back( &( temsp[i][j] ) );
	    }
	    toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<divergenceResults[i][j].printWithJacknife( &jacknivesToSend )<<endl;
	}
      }

      return toReturn.str();
    
}


void SumStatAvgCoa::computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined){
    //currentRow = dataToUse->at(i);
    //cout<<"coord1 "<<dataToUse->at(d).coordinate<<" "<<d<<" "<<dataToUse->at(d) <<endl;
       
    if(!isResolvedDNA(recordToUse->ref) ){
	cerr<<"SumStateAvgCoa:  pairwiseAvgCoa() Problem for record #"<<"\t"<<recordToUse->chr<<" coordinate = "<<recordToUse->coordinate<<" reference = "<<recordToUse->ref<<" is not resolved"<<endl;
	exit(1);
    }
    if(!isResolvedDNA(recordToUse->alt)){ // if one of A,C,G,T
	return ; //next iteration, we need a valid alternative allele otherwise there is no impact on divergence calculation, this speeds it up
    }

    //first one is the ancestral
    //last one is the human reference
    char sampledAllele[ numberOfPopulations]; //array of sampled alleles
    bool cpgForPop    [ numberOfPopulations]; //array of flags to say if the current record is cpg

    //initialize the sampledAllele and cpgForPop

    //double check, already checked in GlacParser 
    if(recordToUse->vectorAlleles->size() != (numberOfPopulations-1)){
	cerr<<"SumStateAvgCoa:  pairwiseAvgCoa() Problem for line "<<recordToUse->chr<<" "<<recordToUse->coordinate<<" wrong number of columns"<<endl;
	exit(1);
    }

    // cout<<"state 2"<<endl;
    // cout<<"coord "<<recordToUse->coordinate<<endl;

    //start at 1 for ancestral
    for(unsigned i=1;i<recordToUse->vectorAlleles->size();i++){
	if(i == 1 ){ //the root can be absent, need to check	      
	    //if the allele count is unknown for both, skip
	    if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
	       recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
		//goto SKIPTONEXTITERATION;
		return ;
	    }
	} 

	if(allowUndefined){
	    if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
	       recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
		sampledAllele[i] = 'N';
		cpgForPop[i]     = false;
		continue ;
	    }
	}

	//plus one for the human allele in sampledAllele and cpgForPop
	sampledAllele[i] =  sampleRandomRefAltAllele(recordToUse->ref,recordToUse->alt,
						     recordToUse->vectorAlleles->at(i).getRefCount(),
						     recordToUse->vectorAlleles->at(i).getAltCount());
	cpgForPop[i] = recordToUse->vectorAlleles->at(i).getIsCpg();
    }

    //storing the human refernce
    sampledAllele[numberOfPopulations-1] = recordToUse->ref;
    cpgForPop[ numberOfPopulations-1]    = recordToUse->vectorAlleles->at(0).getIsCpg();//set the cpg to the ancestral CpG flag

    //debug
    // cout<<"state 2"<<endl;
    // cout<<recordToUse.coordinate<<endl;
    // for(unsigned int ind=0;ind<numberOfPopulations;ind++)
    // 	   cout<<"sampledAllele["<<ind<<"]\t"<<sampledAllele[ind]<<endl;
    // //exit(1);


    //the first is the ancestral, we should never reach that state given the check above
    if(sampledAllele[1] == 'N'){//the ancestral has an unresolved allele, skip
	//goto SKIPTONEXTITERATION;
	//continue;
	return ;
    }
    // cout<<"record "<<alleleRecordsAsString(*currentRow)<<endl;
    //segregatingSites.push_back(*currentRow);//copy constructor




    //for each population, except the root/ancestral at index 0,1
    for(unsigned i=2;i<numberOfPopulations;i++){

	//for each population, except the root/ancestral at index 0,1
	for(unsigned j=2;j<numberOfPopulations;j++){	       
	    //skip when the allele is identical
	    if(i==j)
		continue;
	       	      

	    if(allowUndefined){//if one has undefined allele
		if(sampledAllele[i] == 'N')
		    continue;
		if(sampledAllele[j] == 'N')
		    continue;
	    }

	    computeDiv(sampledAllele[1], //0 is root, 1 is ancestral
		       sampledAllele[i],
		       sampledAllele[j],
		       (cpgForPop[i] || cpgForPop[j]),
		       &(divergenceResults[i][j])  );
	       
	}
    }



    // //for each population, except the root/ancestral at index 0,1
    // for(unsigned i=2;i<numberOfPopulations;i++){

    // 	//for each population, except the root/ancestral at index 0,1
    // 	for(unsigned j=2;j<numberOfPopulations;j++){	       
    // 	    //skip when the allele is identical
    // 	    if(i==j)
    // 		continue;
	       	      

    // 	    if(allowUndefined){//if one has undefined allele
    // 		if(sampledAllele[i] == 'N')
    // 		    continue;
    // 		if(sampledAllele[j] == 'N')
    // 		    continue;
    // 	    }

    // 	    computeDiv(sampledAllele[1], //0 is root, 1 is ancestral
    // 		       sampledAllele[i],
    // 		       sampledAllele[j],
    // 		       (cpgForPop[i] || cpgForPop[j]),
    // 		       &(divergenceResults[i][j])  );
	       
    // 	}
    // }

    //SKIPTONEXTITERATION:
    return;
    
}



void SumStatAvgCoa::computeStat( const   vector < AlleleRecords >  * dataToUse,const vector<string> * popNames,const bool allowUndefined){
    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    divergenceResults = new AvgCoaResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
   }
   populationNames->push_back("ref");

   //vector<AlleleRecords> segregatingSites;

   //while(mp.hasData()){
   //AlleleRecords * currentRow;
   for(unsigned int d=0;d<dataToUse->size();d++){
       computeStatSingle(&(dataToUse->at(d)),allowUndefined);
	  
   }//while the parser has data
   //cout<<"done Sum"<<endl;
   // for(unsigned int i=0;i<numberOfPopulations;i++)
   //     for(unsigned int j=0;j<numberOfPopulations;j++)
   // 	   cout<<i<<" "<<j<<" "<<divergenceResults[i][j]<<endl;
   ///++)
   //   delete(populationNames);

}

void SumStatAvgCoa::read(const string & res){
    //cout<<"READ() SSC"<<endl;
    //return s;
    vector<string> lines = allTokens(res,'\n');

    //DETECTING # OF POP/IND
    populationNames = new vector<string>();
    map<string,int> ind2index;
    unsigned int indIndices=2;//indIndices =0 (root)  indIndices =1 (anc) 
    ind2index["root"]=0;  populationNames->push_back( "root" );
    ind2index["anc"]=1;   populationNames->push_back( "anc" );


    for(unsigned int i=0;i<lines.size();i++){
	if(strBeginsWith(lines[i],"---")) continue; //separator
	if( lines[i].empty() ) continue; //separator
	//cerr<<lines[i]<<endl;
	//vector<string> fields = allTokens(lines[i],'\t');
	string indsIds;	
	for(unsigned int j=0;j<lines[i].size();j++){
	    if( lines[i][j] == '\t') break;
	    indsIds+=lines[i][j];
	}

	string ind1;
	string ind2;
	unsigned int k=0;
	while(k<indsIds.size()){
	    if(indsIds[k] == '-') break;
	    ind1+=indsIds[k++];  
	}
	k++;
	while(k<indsIds.size()){
	    if(indsIds[k] == '\t') break;
	    ind2+=indsIds[k++];	   
	}

	//cout<<ind1<<" "<<ind2<<" "<<src<<endl;
	trimWhiteSpacesBothEnds(&ind1);
	trimWhiteSpacesBothEnds(&ind2);

	if(ind2index.find(ind1) == ind2index.end()){
	    ind2index[ind1]=indIndices;
	    populationNames->push_back(  ind1 );
	    indIndices++;
	}

	if(ind2index.find(ind2) == ind2index.end()){
	    ind2index[ind2]=indIndices;
	    populationNames->push_back(  ind2 );
	    indIndices++;
	}
	

	//cout<<fields.size()<<endl;
    }
    // cerr<<"indIndices "<<indIndices<<endl;
    // for( map<string,int>::const_iterator it=ind2index.begin();it!=ind2index.end();it++){
    // 	cerr<<it->first<<"\t"<<it->second<<endl;
    // }
    numberOfPopulations=indIndices;

    divergenceResults = new AvgCoaResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];


    for(unsigned int i=0;i<lines.size();i++){
	if(strBeginsWith(lines[i],"---")) continue; //separator
	if( lines[i].empty() )            continue; //separator

	//vector<string> fields = allTokens(lines[i],'\t');
	//cout<<"line read "<<lines[i]<<endl;
	string indsIds;	
	for(unsigned int j=0;j<lines[i].size();j++){
	    if( lines[i][j] == '\t') break;
	    indsIds+=lines[i][j];
	}
	
	string ind1;
	string ind2;
	
	unsigned int k=0;
	while(k<indsIds.size()){
	    if(indsIds[k] == '-')		break;	   	    
	    ind1+=indsIds[k++];  
	}
	k++;
	//cerr<<"#l "<<ind1<<"# #"<<ind2<<"# "<<endl;
	while(k<indsIds.size()){
	    if(indsIds[k] == '\t')		break;	    
	    ind2+=indsIds[k++];	   
	}
	k++;
	trimWhiteSpacesBothEnds(&ind1);
	trimWhiteSpacesBothEnds(&ind2);
	//cerr<<"#m "<<ind1<<"# #"<<ind2<<"# "<<endl;
	//fields.pop();//remove first element
	//cout<<(lines[i]+indsIds.size() )<<endl;
	//cout<<"#"<<lines[i].substr(indsIds.size()+1,lines[i].size())<<"#"<<endl;
	istringstream in ( lines[i].substr(indsIds.size()+1,lines[i].size()) );
	// cerr<<endl<<"line  "<<lines[i].substr(indsIds.size()+1,lines[i].size())<<endl;
	// cerr<<"index "<<ind1<<" "<<ind2<<" " <<endl;
	// cerr<<"index "<<ind2index[ind1]<<" "<<ind2index[ind2]<<endl;
	//lines[i].substr(indsIds.size()+1,lines[i].size())<<endl;

	//                    @[src]              [ind1]       -     [ind2]
	//in>>dstatResults[ ind2index[src] ] [ ind2index[ind1] ][ ind2index[ind2] ];
	in>>divergenceResults[ ind2index[ind1] ][ ind2index[ind2] ];
	//cerr<<"#stat "<<f2Results[ ind2index[ind1] ][ ind2index[ind2] ]<<endl;
	//cout<<fields.size()<<endl;
    }

    // istringstream in (strResults);
    //     dstatResults[i][j][k].
    //[numberOfPopulations][numberOfPopulations];
    //storing results
    
   //for(unsigned int i=0;i<popNames->size();i++){
   //	populationNames->push_back(  popNames->at(i) );
   //}
   //populationNames->push_back("ref");
    
   //vector<AlleleRecords> segregatingSites;

   //while(mp.hasData()){


}

string SumStatAvgCoa::print() const {
    // cout<<"print() begin"<<endl;
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){
       for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   
	   toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<divergenceResults[i][j]<<endl;
       }
   }
    // cout<<"print() done"<<endl;
    return toReturn.str();
}

