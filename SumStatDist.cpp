/*
 * SumStatDist
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SumStatDist.h"


SumStatDist::SumStatDist(){

}

SumStatDist::SumStatDist(const SumStatDist & other){
    // cout<<"copy"<<endl;
    numberOfPopulations = other.numberOfPopulations;
    // cout<<"copy "<<numberOfPopulations<<endl;

    populationNames = new vector<string>();

    for(unsigned int i=0;i<other.populationNames->size();i++){ 
	populationNames->push_back(other.populationNames->at(i));
	// cout<<"copyNames "<< other.populationNames->at(i) <<endl;
    }

    divergenceResults = new DistResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new DistResult[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
	    divergenceResults[i][j]=other.divergenceResults[i][j];	    
	    // cout<<"copy "<<divergenceResults[i][j]<<endl;
	}
    }
    // cout<<"copy"<<endl;
}


SumStatDist::SumStatDist(const vector<string> * popNames){

    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    divergenceResults = new DistResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new DistResult[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
    }
    populationNames->push_back("href");
    
}

// SumStatDist::SumStatDist(const SumStatDist & other){

// }

SumStatDist::~SumStatDist(){
    // cout<<"destructor#1"<<endl;
    delete(populationNames);
    //cout<<"destructor#2"<<endl;
    for(unsigned int i=0;i<numberOfPopulations;i++){
	//cout<<"destructor #"<<i<<endl;
	delete [] divergenceResults[i];
    }
    delete [] divergenceResults;
}

DistResult const * const *  SumStatDist::getDistResult() const{
    return divergenceResults;
}

string SumStatDist::printWithBootstraps(const   vector<SumStatDist *> * jackVec) const{ //all boostraps or jacknife
    // string toReturn;
    
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){

	//for each population, except the root/ancestral at index 0,1
	for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   

	    vector< const DistResult * > jacknivesToSend;
	    for(unsigned k=0;k<jackVec->size();k++){	       
		const DistResult * const* temsp= (jackVec->at(k)->getDistResult());
		//constDistResult
		//cout<<temsp[i][j]<<endl;
		//const DistResult * test1  = &temsp[i][j];
		//jacknivesToSend.push_back( &(  (jackVec->at(k)->getDistResult())[i][j]) );
		jacknivesToSend.push_back( &( temsp[i][j] ) );
	    }
	    toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<divergenceResults[i][j].printWithJacknife( &jacknivesToSend )<<endl;
	}
      }

      return toReturn.str();
    
}


void SumStatDist::computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined){
    //currentRow = dataToUse->at(i);
    //cout<<"coord1 "<<dataToUse->at(d).coordinate<<" "<<d<<" "<<dataToUse->at(d) <<endl;
       
    if(!isResolvedDNA(recordToUse->ref) ){
	cerr<<"SumStateDist:  pairwiseDist() Problem for record #"<<"\t"<<recordToUse->chr<<" coordinate = "<<recordToUse->coordinate<<" reference = "<<recordToUse->ref<<" is not resolved"<<endl;
	exit(1);
    }

    //first one is the ancestral
    //last one is the human reference
    char sampledAllele[ numberOfPopulations]; //array of sampled alleles
    bool cpgForPop    [ numberOfPopulations]; //array of flags to say if the current record is cpg

    //initialize the sampledAllele and cpgForPop
    if(!isResolvedDNA(recordToUse->alt)){ // if one of A,C,G,T
	return ; //next iteration, we need a valid alternative allele otherwise there is no impact on divergence calculation, this speeds it up
    }

    //double check, already checked in GlacParser 
    if(recordToUse->vectorAlleles->size() != (numberOfPopulations-1)){
	cerr<<"SumStateDist:  pairwiseDist() Problem for line "<<recordToUse->chr<<" "<<recordToUse->coordinate<<" wrong number of columns"<<endl;
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



void SumStatDist::computeStat( const   vector < AlleleRecords >  * dataToUse,const vector<string> * popNames,const bool allowUndefined){
    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    divergenceResults = new DistResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	divergenceResults[i] = new DistResult[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
   }
   populationNames->push_back("href");

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


string SumStatDist::print() const {
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

