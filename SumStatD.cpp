/*
 * SumStatD
 * Date: Oct-14-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SumStatD.h"


SumStatD::SumStatD(){
    //cout<<"CONSTR1"<<endl;
}

SumStatD::SumStatD(const SumStatD & other){
    // cout<<"copy"<<endl;
    numberOfPopulations = other.numberOfPopulations;
    // cout<<"copy "<<numberOfPopulations<<endl;

    populationNames = new vector<string>();

    for(unsigned int i=0;i<other.populationNames->size();i++){ 
	populationNames->push_back(other.populationNames->at(i));
	// cout<<"copyNames "<< other.populationNames->at(i) <<endl;
    }


    dstatResults = new DstatResult**[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
     	dstatResults[i] = new DstatResult*[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
    	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
    	    dstatResults[i][j] = new DstatResult[numberOfPopulations];	    
    	}
    }

    for(unsigned i=0;i<numberOfPopulations;i++){
    	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
	    for(unsigned k=0;k<numberOfPopulations;k++){	       	    
		dstatResults[i][j][k] = other.dstatResults[i][j][k];	    
	    }
    	}
    }


    // divergenceResults = new AvgCoaResult*[numberOfPopulations];
    // for(unsigned int i=0;i<numberOfPopulations;i++)
    // 	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];

    // for(unsigned i=0;i<numberOfPopulations;i++){
    // 	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
    // 	    divergenceResults[i][j]=other.divergenceResults[i][j];	    
    // 	    // cout<<"copy "<<divergenceResults[i][j]<<endl;
    // 	}
    // }
    // cout<<"copy"<<endl;

}


SumStatD::SumStatD(const vector<string> * popNames){
    
    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    // divergenceResults = new AvgCoaResult*[numberOfPopulations];
    // for(unsigned int i=0;i<numberOfPopulations;i++)
    // 	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];
    // //[numberOfPopulations][numberOfPopulations];


    dstatResults = new DstatResult**[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
     	dstatResults[i] = new DstatResult*[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
    	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
    	    dstatResults[i][j] = new DstatResult[numberOfPopulations];	    
    	}
    }


    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
    }
    populationNames->push_back("ref");
    
}

// SumStatD::SumStatD(const SumStatD & other){

// }

SumStatD::~SumStatD(){
  // cout<<"destructor#1"<<endl;
    delete(populationNames);
    //cout<<"destructor#2"<<endl;
    //cout<<"destructor #"<<i<<endl;

    for(unsigned int i=0;i<numberOfPopulations;i++){
	for(unsigned int j=0;j<numberOfPopulations;j++){
	    delete [] dstatResults[i][j];
	}
    }

    for(unsigned int i=0;i<numberOfPopulations;i++){
	delete [] dstatResults[i];
    }

    delete [] dstatResults;
}


DstatResult const * const * const * SumStatD::getDstatResult() const{
    return dstatResults;
}

string SumStatD::printWithBootstraps(const   vector<SumStatD *> * jackVec, const string & dnaDistMode) const{ //all boostraps or jacknife
    // string toReturn;
    
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){

	//for each population, except the root/ancestral at index 0,1
	for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   

	   for(unsigned k=2;k<numberOfPopulations;k++){	       
	       //skip when the population is identical
	       if(i==k)
		   continue;	   
	       if(j==k)
		   continue;	   

	       vector< const DstatResult * > jacknivesToSend;
	       for(unsigned l=0;l<jackVec->size();l++){	       
		   const DstatResult * const*  const* temsp = (jackVec->at(l)->getDstatResult());
		   //constAvgCoaResult
		   //cout<<temsp[i][j]<<endl;
		   //const AvgCoaResult * test1  = &temsp[i][j];
		   //jacknivesToSend.push_back( &(  (jackVec->at(k)->getAvgCoaResult())[i][j]) );
		   jacknivesToSend.push_back( &( temsp[i][j][k] ) );
	       }
	       //toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"\t"<<divergenceResults[i][j][k].printWithJacknife( &jacknivesToSend )<<endl;
               toReturn<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"@"<< populationNames->at(i)  <<"\t"<<dstatResults[i][j][k].printWithJacknife( &jacknivesToSend )<<endl;

	   }//k
	}//j
    }//i

    return toReturn.str();

}


void SumStatD::computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined){
    //currentRow = dataToUse->at(i);
    //cout<<"coord1 "<<dataToUse->at(d).coordinate<<" "<<d<<" "<<dataToUse->at(d) <<endl;
    //cout<<"coord1 "<<recordToUse->coordinate<<endl;
       
    if(!isResolvedDNA(recordToUse->ref) ){
	cerr<<"SumStatD  computeStatSingle() Problem for record #"<<"\t"<<recordToUse->chr<<" coordinate = "<<recordToUse->coordinate<<" reference = "<<recordToUse->ref<<" is not resolved"<<endl;
	exit(1);
    }
    //initialize the sampledAllele and cpgForPop
    if(!isResolvedDNA(recordToUse->alt)){ // if one of A,C,G,T
	return ; //next iteration, we need a valid alternative allele otherwise there is no impact on divergence calculation, this speeds it up
    }

    //first one is the ancestral
    //last one is the human reference
    char sampledAllele[ numberOfPopulations]; //array of sampled alleles
    bool cpgForPop    [ numberOfPopulations]; //array of flags to say if the current record is cpg


    //double check, already checked in GlacParser 
    if(recordToUse->vectorAlleles->size() != (numberOfPopulations-1)){
	cerr<<"SumStatD.cpp  computeStatSingle() Problem for line "<<recordToUse->chr<<" "<<recordToUse->coordinate<<" wrong number of columns"<<endl;
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
    // 	cout<<"sampledAllele["<<ind<<"]\t"<<sampledAllele[ind]<<endl;
    // //exit(1);


    //the first is the ancestral, we should never reach that state given the check above
    if(sampledAllele[1] == 'N'){//the ancestral has an unresolved allele, skip
	//goto SKIPTONEXTITERATION;
	//continue;
	return ;
    }
    // cout<<"record "<<alleleRecordsAsString(*currentRow)<<endl;
    //segregatingSites.push_back(*currentRow);//copy constructor

    //cout<<"coord1 "<<recordToUse->coordinate<<endl;

    //for each population, except the root/ancestral at index 0,1
    for(unsigned i=2;i<numberOfPopulations;i++){

	//for each population, except the root/ancestral at index 0,1
	for(unsigned j=2;j<numberOfPopulations;j++){	       
	    //skip when the allele is identical
	    if(i==j)
		continue;

	    for(unsigned k=2;k<numberOfPopulations;k++){	       
		//skip when the allele is identical
		if(i==k)
		    continue;
		if(j==k)
		    continue;
	       	      

		if(allowUndefined){//if one has undefined allele
		    if(sampledAllele[i] == 'N')
			continue;
		    if(sampledAllele[j] == 'N')
			continue;
		    if(sampledAllele[k] == 'N')
			continue;
		}
		//cout<<"seg "<<recordToUse->coordinate<<"\t"<<i<<"\t"<<j<<"\t"<<k<<endl;
		//bool dstval = computeDstat(sampledAllele[0], //root
		computeDstat(sampledAllele[1], //root
			     sampledAllele[i], //derived
			     sampledAllele[j], //ind 1
			     sampledAllele[k], //ind 2
			     (cpgForPop[j] || cpgForPop[k]), //only look at j and k for CpG
			     &(dstatResults[i][j][k]) );
		//cout<<dstatResults[i][j][k]<<endl;
		// computeDiv(sampledAllele[1], //0 is root, 1 is ancestral
		// 	   sampledAllele[i],
		// 	   sampledAllele[j],
		// 	   (cpgForPop[i] || cpgForPop[j]),
		// 	   &(divergenceResults[i][j])  );
	    }//k
	}//j
    }//i

    
    
    //SKIPTONEXTITERATION:
    return;

}



void SumStatD::computeStat( const   vector < AlleleRecords >  * dataToUse,const vector<string> * popNames,const bool allowUndefined){
    // cout<<"computeStat() "<<endl;
   numberOfPopulations=popNames->size()+1;//+1 for the human reference
   dstatResults = new DstatResult**[numberOfPopulations];

    for(unsigned int i=0;i<numberOfPopulations;i++)
	dstatResults[i] = new DstatResult*[numberOfPopulations];

    for(unsigned int i=0;i<numberOfPopulations;i++)
	for(unsigned int j=0;j<numberOfPopulations;j++)
	    dstatResults[i][j] = new DstatResult[numberOfPopulations];

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


string SumStatD::print() const {
    // cout<<"SumStatD print() begin"<<endl;
    // exit(1);
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){
       for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   
	    for(unsigned k=2;k<numberOfPopulations;k++){	       
		//skip when the allele is identical
		if(i==k)
		    continue;
		if(j==k)
		    continue;
		// cout<<"print() "<<i<<" "<<j<<" "<<k<<endl;
		// cout<<"1#"<<populationNames->at(j)<<endl;
		// cout<<"2#"<<populationNames->at(k)<<endl;
		// cout<<"3#"<<populationNames->at(i)  <<endl;
		// cout<<"4#"<<dstatResults[i][j][k]<<endl;
		toReturn<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"@"<< populationNames->at(i)  <<"\t"<<dstatResults[i][j][k]<<endl;
		//toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<divergenceResults[i][j]<<endl;
	    }
       }
   }
    // cout<<"print() done"<<endl;
    return toReturn.str();
}

