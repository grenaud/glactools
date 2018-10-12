/*
 * SumStatF3
 * Date: Oct-14-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SumStatF3.h"


SumStatF3::SumStatF3(){
    //cout<<"CONSTR1"<<endl;
}

SumStatF3::SumStatF3(const SumStatF3 & other){
    // cout<<"copy"<<endl;
    numberOfPopulations = other.numberOfPopulations;
    // cout<<"copy "<<numberOfPopulations<<endl;

    populationNames = new vector<string>();

    for(unsigned int i=0;i<other.populationNames->size();i++){ 
	populationNames->push_back(other.populationNames->at(i));
	// cout<<"copyNames "<< other.populationNames->at(i) <<endl;
    }


    f3Results = new F3Result**[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
     	f3Results[i] = new F3Result*[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
    	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
    	    f3Results[i][j] = new F3Result[numberOfPopulations];	    
    	}
    }

    for(unsigned i=0;i<numberOfPopulations;i++){
    	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
	    for(unsigned k=0;k<numberOfPopulations;k++){	       	    
		f3Results[i][j][k] = other.f3Results[i][j][k];	    
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


SumStatF3::SumStatF3(const vector<string> * popNames){
    
    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    // divergenceResults = new AvgCoaResult*[numberOfPopulations];
    // for(unsigned int i=0;i<numberOfPopulations;i++)
    // 	divergenceResults[i] = new AvgCoaResult[numberOfPopulations];
    // //[numberOfPopulations][numberOfPopulations];



    f3Results = new F3Result**[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
     	f3Results[i] = new F3Result*[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
    	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
    	    f3Results[i][j] = new F3Result[numberOfPopulations];	    
    	}
    }


    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
    }
    populationNames->push_back("href");
    
}

// SumStatF3::SumStatF3(const SumStatF3 & other){

// }

SumStatF3::~SumStatF3(){
  // cout<<"destructor#1"<<endl;
    delete(populationNames);
    //cout<<"destructor#2"<<endl;
    //cout<<"destructor #"<<i<<endl;

    for(unsigned int i=0;i<numberOfPopulations;i++){
	for(unsigned int j=0;j<numberOfPopulations;j++){
	    delete [] f3Results[i][j];
	}
    }

    for(unsigned int i=0;i<numberOfPopulations;i++){
	delete [] f3Results[i];
    }

    delete [] f3Results;
}


F3Result const * const * const * SumStatF3::getF3Result() const{
    return f3Results;
}

string SumStatF3::printWithBootstraps(const   vector<SumStatF3 *> * jackVec, const string & dnaDistMode) const{ //all boostraps or jacknife
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

	       vector< const F3Result * > jacknivesToSend;
	       for(unsigned l=0;l<jackVec->size();l++){	       
		   const F3Result * const*  const* temsp = (jackVec->at(l)->getF3Result());
		   //constAvgCoaResult
		   //cout<<temsp[i][j]<<endl;
		   //const AvgCoaResult * test1  = &temsp[i][j];
		   //jacknivesToSend.push_back( &(  (jackVec->at(k)->getAvgCoaResult())[i][j]) );
		   jacknivesToSend.push_back( &( temsp[i][j][k] ) );
	       }
	       //toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"\t"<<divergenceResults[i][j][k].printWithJacknife( &jacknivesToSend )<<endl;
               toReturn<<populationNames->at(j)<<"-"<<populationNames->at(k)<<"@"<< populationNames->at(i)  <<"\t"<<f3Results[i][j][k].printWithJacknife( &jacknivesToSend )<<endl;

	   }//k
	}//j
    }//i

    return toReturn.str();

}


void SumStatF3::computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined){
    //currentRow = dataToUse->at(i);
    //cout<<"coord1 "<<dataToUse->at(d).coordinate<<" "<<d<<" "<<dataToUse->at(d) <<endl;
    //cerr<<"coord1 "<<recordToUse->coordinate<<endl;
       
    if(!isResolvedDNA(recordToUse->ref) ){
	cerr<<"SumStatF3  computeStatSingle() Problem for record #"<<"\t"<<recordToUse->chr<<" coordinate = "<<recordToUse->coordinate<<" reference = "<<recordToUse->ref<<" is not resolved"<<endl;
	exit(1);
    }

    if(!isResolvedDNA(recordToUse->alt)){ // if one of A,C,G,T
	return ; //next iteration, we need a valid alternative allele otherwise there is no impact on divergence calculation, this speeds it up
    }

    //first one is the ancestral
    //todo remove
    char sampledAllele[   numberOfPopulations]; //array of sampled alleles
    double freqAllele [   numberOfPopulations]; //array of sampled alleles    
    int  refCount [   numberOfPopulations]; 
    int  altCount [   numberOfPopulations]; 

    bool cpgForPop    [   numberOfPopulations]; //array of flags to say if the current record is cpg
    bool undefined    [   numberOfPopulations]; //array of flags to say if the current record is cpg

    //initialize the sampledAllele and cpgForPop

    //double check, already checked in GlacParser 
    if(recordToUse->vectorAlleles->size() != (numberOfPopulations-1)){
	cerr<<"SumStatF3.cpp  computeStatSingle() Problem for line "<<recordToUse->chr<<" "<<recordToUse->coordinate<<" wrong number of columns"<<endl;
	exit(1);
    }

    // cout<<"state 2"<<endl;
    //cerr<<"coord "<<recordToUse->coordinate<<endl;
    if( recordToUse->vectorAlleles->at(1).getRefCount() == 0 &&
	recordToUse->vectorAlleles->at(1).getAltCount() == 0){
	return ;
    }

    char ancBase =sampleRandomRefAltAllele(recordToUse->ref,
					   recordToUse->alt,
					   recordToUse->vectorAlleles->at(1).getRefCount(),
					   recordToUse->vectorAlleles->at(1).getAltCount());
    
    bool refIsAnc= (ancBase == recordToUse->ref);
    
    //start at 1 for ancestral
    // unsigned int derAllele=0;
    // unsigned int ancAllele=0;
    unsigned int refAllele=0;
    unsigned int altAllele=0;
    
    bool isSitePotentialTransition ;
    bool isSitePotentialDamage     ;
    double m;
    double sumF=0.0;
    for(unsigned i=1;i<recordToUse->vectorAlleles->size();i++){
	// if(i == 1 ){ //the root can be absent, need to check	      
	//     //if the allele count is unknown for both, skip
	//     if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
	//        recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
	// 	//goto SKIPTONEXTITERATION;
	// 	return ;
	//     }
	// } 
	
	if(!allowUndefined){
	    if(recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&
	       recordToUse->vectorAlleles->at(i).getAltCount() ==  0){
		// 	//sampledAllele[i] = 'N';		
		// 	cpgForPop[i]     = false;
		// 	continue ;
		goto SKIPTONEXTITERATION;
	    }
	}

	refAllele+=recordToUse->vectorAlleles->at(i).getRefCount();
	altAllele+=recordToUse->vectorAlleles->at(i).getAltCount();
	//taken from CountData treemix
	double sumAllBases = double( recordToUse->vectorAlleles->at(i).getRefCount() + recordToUse->vectorAlleles->at(i).getAltCount() );
	double freqRef = double( recordToUse->vectorAlleles->at(i).getRefCount() ) / sumAllBases;
	freqAllele[i] =  freqRef;
	sumF+=freqRef;
	// average_nInds[ j ] += sumAllBases/2.0;
	
	// double tmp2 = (double) recordToUse->vectorAlleles->at(i).getAltCount() / (sumAllBases - 1.0);
	// double tmphzy = 2* freqRef * tmp2;
	// if ( sumAllBases < 2){
	//     tmphzy = 2*freqRef*(1-freqRef);
	// }
	// average_hzy[j] += tmphzy; //2*f*(1-f);
	// id2nsnp++;


#ifdef WEIRDF
	if(refIsAnc){//alt is derived
	    derAllele+=recordToUse->vectorAlleles->at(i).getAltCount();
	    ancAllele+=recordToUse->vectorAlleles->at(i).getRefCount();
	    freqAllele[i] = double( recordToUse->vectorAlleles->at(i).getAltCount() ) / double( recordToUse->vectorAlleles->at(i).getRefCount() + recordToUse->vectorAlleles->at(i).getAltCount() );
	}else{//the ref is derived, alt is ancestral	    
	    derAllele+=recordToUse->vectorAlleles->at(i).getRefCount();
	    ancAllele+=recordToUse->vectorAlleles->at(i).getAltCount();
	    freqAllele[i] = double( recordToUse->vectorAlleles->at(i).getRefCount() ) / double( recordToUse->vectorAlleles->at(i).getRefCount() + recordToUse->vectorAlleles->at(i).getAltCount() );	    
	}
#endif

	refCount [ i ] = recordToUse->vectorAlleles->at(i).getRefCount();
	altCount [ i ] = recordToUse->vectorAlleles->at(i).getAltCount();
	//cerr<<i<<" f="<< freqAllele[i]<<endl;
	//plus one for the human allele in sampledAllele and cpgForPop
	// sampledAllele[i] =  sampleRandomRefAltAllele(recordToUse->ref,recordToUse->alt,
	// 					     recordToUse->vectorAlleles->at(i).getRefCount(),
	// 					     recordToUse->vectorAlleles->at(i).getAltCount());
	cpgForPop[i] = recordToUse->vectorAlleles->at(i).getIsCpg();
	undefined[i]   = (recordToUse->vectorAlleles->at(i).getRefCount() ==  0 &&        recordToUse->vectorAlleles->at(i).getAltCount() ==  0);
	//cerr<<i<<" "<<recordToUse->vectorAlleles->at(i).getRefCount()<<" "<< recordToUse->vectorAlleles->at(i).getAltCount() <<" "<<freqAllele[i]<<endl;
    }
    //we do not add the refAllele into the count


    //storing the human refernce
#ifdef WEIRDF
    if(refIsAnc){//alt is derived
    	freqAllele[numberOfPopulations-1] = double( 0 ) / double( 1 + 0 );  //100% re
    }else{//ref is der, alt is anc
    	freqAllele[numberOfPopulations-1] = double( 1 ) / double( 1 + 0 );  // 0% ref
    }
#endif



    freqAllele[numberOfPopulations-1] = double( 1 ) / double( 1 + 0 );  // 0% ref    

    refCount [ numberOfPopulations-1 ] = 1;
    altCount [ numberOfPopulations-1 ] = 0;

    // sampledAllele[numberOfPopulations-1] = recordToUse->ref;
    cpgForPop[ numberOfPopulations-1]    = recordToUse->vectorAlleles->at(0).getIsCpg();//set the cpg to the ancestral CpG flag
    undefined[ numberOfPopulations-1 ]   = false;//we always have the ref

    //m = double(refAllele)/double(refAllele+altAllele);
    m = double(sumF)/double(numberOfPopulations-2);//we do not count the ref+(either root or anc)

    //cerr<<"m "<<m<<" "<<numberOfPopulations<<endl;
    // if(randomBool()){//flip the allele frequencies at each 0.5 site
	
    // 	for(unsigned i=2;i<numberOfPopulations;i++)
    // 	    freqAllele [ i ] = 1.0 - freqAllele [ i ]; 
    // }
    
    //if(true)//TODO to re-enable??
    for(unsigned i=1;i<(recordToUse->vectorAlleles->size()+1);i++){//go over by one to compute the reference
	//cerr<<i<<" f="<< freqAllele[i]<<" m="<< m<<endl;
	freqAllele[i] =   freqAllele[i] - m;
	//cerr<<i<<" f="<< freqAllele[i]<<"\t"<<undefined[i]<<endl;
	// average_nInds[i] = average_nInds[i]/ id2nsnp[i];
	// average_hzy[i]   = average_hzy[i]/ id2nsnp[i];
    }
    
    
    isSitePotentialTransition = isPotentialTransition(  recordToUse->ref, recordToUse->alt );
    isSitePotentialDamage     = isSitePotentialTransition;//since for f3 we do not care about the ancestral, we cannot tell anything. 
    

    //update base counters
    


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
	       	      

		if(!allowUndefined){//if one has undefined allele
		    if(undefined[i] )
			continue;
		    if(undefined[j] )
			continue;
		    if(undefined[k] )
			continue;
		}
		//cerr<<allowUndefined<<" seg "<<recordToUse->coordinate<<"\t"<<i<<"\t"<<j<<"\t"<<k<<endl;
		//bool dstval = computeF3(sampledAllele[0], //root
		//cerr<<"precounter "<<i<<" "<<j<<" "<<k<<endl;
		addCountersInd(     refCount [ i ] ,
				    altCount [ i ] ,
				    (cpgForPop[i] || cpgForPop[j] || cpgForPop[k]), //only look at j and k for CpG
				    isSitePotentialTransition,
				    isSitePotentialDamage,
				    &(f3Results[i][j][k]) );
		//cerr<<"pstcounter "<<i<<" "<<j<<" "<<k<<endl;

		computeF3triple(//freqAllele[1], //root
			  // (freq_condition-freq_ind1)*(freq_condition-freq_ind2);
			  freqAllele[i], //condition or target
			  freqAllele[j], //ind 1
			  freqAllele[k], //ind 2
			  (cpgForPop[i] || cpgForPop[j] || cpgForPop[k]), //only look at j and k for CpG
			  isSitePotentialTransition,
			  isSitePotentialDamage,
			  &(f3Results[i][j][k]) );
	    }//k
	}//j
    }//i

    
    
    SKIPTONEXTITERATION:
    return;

}



void SumStatF3::computeStat( const   vector < AlleleRecords >  * dataToUse,const vector<string> * popNames,const bool allowUndefined){
    // cout<<"computeStat() "<<endl;
   numberOfPopulations=popNames->size()+1;//+1 for the human reference
   f3Results = new F3Result**[numberOfPopulations];

    for(unsigned int i=0;i<numberOfPopulations;i++)
	f3Results[i] = new F3Result*[numberOfPopulations];

    for(unsigned int i=0;i<numberOfPopulations;i++)
	for(unsigned int j=0;j<numberOfPopulations;j++)
	    f3Results[i][j] = new F3Result[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
   }
   populationNames->push_back("ref");

    // average_nInds = vector<double> (numberOfPopulations,0.0);
    // average_hzy   = vector<double> (numberOfPopulations,0.0);
    // id2nsnp       = vector<double> (numberOfPopulations,0.0);

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


string SumStatF3::print() const {
    // cout<<"SumStatF3 print() begin"<<endl;
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
		// cout<<"4#"<<f3Results[i][j][k]<<endl;
		toReturn<<populationNames->at(i)<<" ; "<<
		    populationNames->at(j)<<" , "<<populationNames->at(k)  <<"\t"<<f3Results[i][j][k]<<endl;
		//toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<divergenceResults[i][j]<<endl;
	    }
       }
   }
    // cout<<"print() done"<<endl;
    return toReturn.str();
}

