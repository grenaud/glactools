/*
 * SumStatPi
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SumStatPi.h"


SumStatPi::SumStatPi(){

}

SumStatPi::SumStatPi(const SumStatPi & other){
    // cout<<"copy"<<endl;
    numberOfPopulations = other.numberOfPopulations;
    // cout<<"copy "<<numberOfPopulations<<endl;

    populationNames = new vector<string>();

    for(unsigned int i=0;i<other.populationNames->size();i++){ 
	populationNames->push_back(other.populationNames->at(i));
	// cout<<"copyNames "<< other.populationNames->at(i) <<endl;
    }

    pianceResults = new PiResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	pianceResults[i] = new PiResult[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
	    pianceResults[i][j]=other.pianceResults[i][j];	    
	    // cout<<"copy "<<pianceResults[i][j]<<endl;
	}
    }
    // cout<<"copy"<<endl;
}


SumStatPi::SumStatPi(const vector<string> * popNames){
    //cerr<<"const "<<"SumStatPi"<<endl;
    numberOfPopulations=popNames->size()+1;//+1 for the reference +2 for root anc

    // cerr<<"numberOfPopulations1 "<<numberOfPopulations<<endl;
    // for(unsigned int i=0;i<popNames->size();i++){
    // 	cerr<<i<<" "<<popNames->at(i)<<endl;
    // }

    pianceResults = new PiResult*[numberOfPopulations];
    //cerr<<"numberOfPopulations2 "<<numberOfPopulations<<endl;
    for(unsigned int i=0;i<numberOfPopulations;i++){
	//cerr<<"numberOfPopulationsi "<<i<<" "<<numberOfPopulations<<endl;
	pianceResults[i] = new PiResult[numberOfPopulations];
    }
    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    // populationNames->push_back("root");
    // populationNames->push_back("anc");

    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
    }
    populationNames->push_back("ref");

    // for(unsigned int i=0;i<populationNames->size();i++){
    // 	cerr<<"popNames["<<i<<"] = "<<populationNames->at(i) <<endl;
    // }
    // exit(1);
}

// SumStatPi::SumStatPi(const SumStatPi & other){

// }

SumStatPi::~SumStatPi(){
    // cout<<"destructor#1"<<endl;
    delete(populationNames);
    //cout<<"destructor#2"<<endl;
    for(unsigned int i=0;i<numberOfPopulations;i++){
	//cout<<"destructor #"<<i<<endl;
	delete [] piResults[i];
    }
    delete [] piResults;
}

PiResult const * const *  SumStatPi::getPiResult() const{
    return pianceResults;
}

string SumStatPi::printWithBootstraps(const   vector<SumStatPi *> * jackVec) const{ //all boostraps or jacknife
    // string toReturn;
    
    stringstream toReturn;
    for(unsigned i=0;i<numberOfPopulations;i++){

	//for each population, even the root/ancestral at index 0,1
	for(unsigned j=0;j<numberOfPopulations;j++){	       
	    //skip when the population is identical
	    if(i==j)
		continue;	   

	    vector< const PiResult * > jacknivesToSend;
	    for(unsigned k=0;k<jackVec->size();k++){	       
		const PiResult * const* temsp= (jackVec->at(k)->getPiResult());
		//constPiResult
		//cout<<temsp[i][j]<<endl;
		//const PiResult * test1  = &temsp[i][j];
		//jacknivesToSend.push_back( &(  (jackVec->at(k)->getPiResult())[i][j]) );
		jacknivesToSend.push_back( &( temsp[i][j] ) );
	    }
	    toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<pianceResults[i][j].printWithJacknife( &jacknivesToSend  )<<endl;
	}
    }
    
    return toReturn.str();    
}


void SumStatPi::computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined=true){
    //currentRow = dataToUse->at(i);
    //cout<<"coord1 "<<dataToUse->at(d).coordinate<<" "<<d<<" "<<dataToUse->at(d) <<endl;
    //cerr<<"coord1 "<<recordToUse->coordinate <<" "<<allowUndefined<<endl;
    
    if(!isResolvedDNA(recordToUse->ref) ){
	cerr<<"SumStatePi:  pairwisePi() Problem for record #"<<"\t"<<recordToUse->chr<<" coordinate = "<<recordToUse->coordinate<<" reference = "<<recordToUse->ref<<" is not resolved"<<endl;
	exit(1);
    }

    bool notASegSite=(recordToUse->alt == 'N');
    int refCountNotRootAnc=0;
    int altCountNotRootAnc=0;

    //first one is the ancestral
    //last one is the human reference
    char sampledAllele[ numberOfPopulations-1]; //array of sampled alleles
    bool cpgForPop    [ numberOfPopulations-1]; //array of flags to say if the current record is cpg
    
    if(notASegSite){
	goto SKIPTONEXTITERATION;
    }
    //initialize the sampledAllele and cpgForPop
    // if(!isResolvedDNA(recordToUse->alt)){ // if one of A,C,G,T
    // 	return ; //next iteration, we need a valid alternative allele otherwise there is no impact on piance calculation, this speeds it up
    // }

    //double check, already checked in GlacParser 
    if(recordToUse->vectorAlleles->size() != (numberOfPopulations-1)){
	cerr<<"SumStatePi:  pairwisePi() Problem for line "<<recordToUse->chr<<" "<<recordToUse->coordinate<<" wrong number of columns"<<recordToUse->vectorAlleles->size()<<" "<<(numberOfPopulations)<<endl;
	exit(1);
    }

    // cerr<<"state 2"<<endl;
    // cerr<<"coord "<<recordToUse->coordinate<<" "<<notASegSite<<endl;

    //start at 1 for ancestral
    if(notASegSite){
	// for(unsigned i=0;i<recordToUse->vectorAlleles->size();i++){
	//     if(allowUndefined){
	// 	if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
	// 	   recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
	// 	    sampledAllele[i] = 'N';
	// 	    cpgForPop[i]     = false;
	// 	    continue ;
	// 	}
	//     }
	//     cpgForPop[i] = recordToUse->vectorAlleles->at(i).getIsCpg();
	// }
    
    }else{//is seg site
	for(unsigned int i=2;i<recordToUse->vectorAlleles->size();i++){
	    // if(i == 1 ){ //the root can be absent, need to check	      
	    //     //if the allele count is unknown for both, skip
	    //     if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
	    //        recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
	    // 	//goto SKIPTONEXTITERATION;
	    // 	return ;
	    //     }
	    // } 
	    
	    if(allowUndefined){
		if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
		   recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
		    sampledAllele[i] = 'N';
		    cpgForPop[i]     = false;
		    continue ;
		}
	    }

	    refCountNotRootAnc+=recordToUse->vectorAlleles->at(i).getRefCount();
	    altCountNotRootAnc+=recordToUse->vectorAlleles->at(i).getAltCount();

	    //plus one for the human allele in sampledAllele and cpgForPop
	    sampledAllele[i] =  sampleRandomRefAltAllele(recordToUse->ref,
							 recordToUse->alt,
							 recordToUse->vectorAlleles->at(i).getRefCount(),
							 recordToUse->vectorAlleles->at(i).getAltCount());
	    cpgForPop[i] = recordToUse->vectorAlleles->at(i).getIsCpg();
	}
    }

    if(	(refCountNotRootAnc == 0) || (altCountNotRootAnc == 0) ){
	goto SKIPTONEXTITERATION;
    }

    // //storing the human refernce
    // sampledAllele[numberOfPopulations-1] = recordToUse->ref;
    // cpgForPop[    numberOfPopulations-1] = recordToUse->vectorAlleles->at(0).getIsCpg();//set the cpg to the ancestral CpG flag
    
    
    //debug
    //    cout<<"state 2"<<endl;
    // cout<<recordToUse->coordinate<<" ";
    // for(unsigned int ind=0;ind<numberOfPopulations;ind++)
    //  	cout<<sampledAllele[ind];
    // cout<<endl;
    // for(unsigned int ind=0;ind<numberOfPopulations;ind++)
    // 	cout<<"sampledAllele["<<ind<<"]\t"<<sampledAllele[ind]<<endl;
    // //exit(1);
    
    
    //the first is the ancestral, we should never reach that state given the check above
    // if(sampledAllele[1] == 'N'){//the ancestral has an unresolved allele, skip
    // 	//goto SKIPTONEXTITERATION;
    // 	//continue;
    // 	return ;
    // }
    // cout<<"record "<<alleleRecordsAsString(*currentRow)<<endl;
    //segregatingSites.push_back(*currentRow);//copy constructor
    



    //for each population, including the root/ancestral at index 0,1
    if(notASegSite){
	
	// int allePairIndex = baseResolved2int(recordToUse->ref)*5;
	// //cerr<<recordToUse->ref<<" "<<allePairIndex<<endl;
	// for(unsigned i=0;i<numberOfPopulations;i++){
	    
	//     //for each population, including the root/ancestral at index 0,1
	//     for(unsigned j=0;j<numberOfPopulations;j++){	       
	// 	//skip when the allele is identical
	// 	if(i==j)
	// 	    continue;
		
	// 	//defined in ComputePi_core
	// 	addIdenticalSite(sampledAllele[1], //0 is root, 1 is ancestral
	// 			 recordToUse->ref,
	// 			 (cpgForPop[i] || cpgForPop[j]),
	// 			 &(pianceResults[i][j]) ,
	// 			 allePairIndex);
		
		
	//     }
	// }
	
    }else{
	for(unsigned i=2;i<numberOfPopulations;i++){
	    
	    //for each population, excluding the root/ancestral at index 0,1
	    for(unsigned j=2;j<numberOfPopulations;j++){	       
		//skip when the allele is identical
		if(i==j)
		    continue;
		
		
		//if(allowUndefined){//if one has undefined allele
		if(sampledAllele[i] == 'N')
		    continue;
		if(sampledAllele[j] == 'N')
		    continue;
		//}
		
		//cout<<"1pop["<<i<<"]:"<<populationNames->at(i)<<"\tpop["<<j<<"]:"<<populationNames->at(j)<<"\t"<<pianceResults[i][j].printAll()<<endl;

		//defined in ComputePi_core
		computePi(//sampledAllele[1], //0 is root, 1 is ancestral
			    sampledAllele[i],
			    sampledAllele[j],
			    (cpgForPop[i] || cpgForPop[j]),
			    &(piResults[i][j])  );

		//cout<<"2pop["<<i<<"]:"<<populationNames->at(i)<<"\tpop["<<j<<"]:"<<populationNames->at(j)<<"\t"<<distanceResults[i][j].printAll()<<endl;
	    
	    }
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
    // 		       &(distanceResults[i][j])  );
	       
    // 	}
    // }

    SKIPTONEXTITERATION:
    return;
    
}



void SumStatDist::computeStat( const   vector < AlleleRecords >  * dataToUse,const vector<string> * popNames,const bool allowUndefined=true){
    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    //cerr<<"numberOfPopulations "<<numberOfPopulations<<endl;
    distanceResults = new DistResult*[numberOfPopulations-2];
    for(unsigned int i=0;i<(numberOfPopulations-2);i++)
	distanceResults[i] = new DistResult[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
	//cerr<<"popNames["<<i<<"] = "<<populationNames->at(i)<<endl;
    }
    //populationNames->push_back("ref");
    // cerr<<"popNames["<<popNames->size()<<"] = "<<populationNames->at(popNames->size())<<endl;
    // exit(1);
    //vector<AlleleRecords> segregatingSites;

   //while(mp.hasData()){
   //AlleleRecords * currentRow;
   for(unsigned int d=0;d<dataToUse->size();d++){
       // cerr<<"record#d "<<dataToUse->at(d)<<endl;
       computeStatSingle( &(dataToUse->at(d)) , allowUndefined);
	  
   }//while the parser has data
   //cout<<"done Sum"<<endl;
   // for(unsigned int i=0;i<numberOfPopulations;i++)
   //     for(unsigned int j=0;j<numberOfPopulations;j++)
   // 	   cout<<i<<" "<<j<<" "<<distanceResults[i][j]<<endl;
   ///++)
   //   delete(populationNames);

}


string SumStatDist::print() const {
    // cout<<"print() begin"<<endl;
    stringstream toReturn;
    for(unsigned i=0;i<numberOfPopulations;i++){
       for(unsigned j=0;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   
	   toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<distanceResults[i][j]<<endl;
       }
   }
    // cout<<"print() done"<<endl;
    return toReturn.str();
}

