/*
 * SumStatFst
 * Date: Feb-26-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SumStatFst.h"


SumStatFst::SumStatFst(){

}

SumStatFst::SumStatFst(const SumStatFst & other){
    // cout<<"copy"<<endl;
    numberOfPopulations = other.numberOfPopulations;
    // cout<<"copy "<<numberOfPopulations<<endl;

    populationNames = new vector<string>();

    for(unsigned int i=0;i<other.populationNames->size();i++){ 
	populationNames->push_back(other.populationNames->at(i));
	// cout<<"copyNames "<< other.populationNames->at(i) <<endl;
    }

    fstResults = new FstResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	fstResults[i] = new FstResult[numberOfPopulations];

    for(unsigned i=0;i<numberOfPopulations;i++){
	for(unsigned j=0;j<numberOfPopulations;j++){	       	    
	    fstResults[i][j]=other.fstResults[i][j];	    
	    // cout<<"copy "<<fstResults[i][j]<<endl;
	}
    }
    // cout<<"copy"<<endl;
}


SumStatFst::SumStatFst(const vector<string> * popNames){

    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    fstResults = new FstResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	fstResults[i] = new FstResult[numberOfPopulations];

    //[numberOfPopulations][numberOfPopulations];

    populationNames = new vector<string>();
    for(unsigned int i=0;i<popNames->size();i++){
	populationNames->push_back(  popNames->at(i) );
    }
    populationNames->push_back("ref");
    
}

// SumStatFst::SumStatFst(const SumStatFst & other){

// }

SumStatFst::~SumStatFst(){
    // cout<<"destructor#1"<<endl;
    delete(populationNames);
    //cout<<"destructor#2"<<endl;
    for(unsigned int i=0;i<numberOfPopulations;i++){
	//cout<<"destructor #"<<i<<endl;
	delete [] fstResults[i];
    }
    delete [] fstResults;
}

FstResult const * const *  SumStatFst::getFstResult() const{
    return fstResults;
}

string SumStatFst::printWithBootstraps(const   vector<SumStatFst *> * jackVec,  const string & dnaDistMode) const{ //all boostraps or jacknife
    // string toReturn;
    
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){

	//for each population, except the root/ancestral at index 0,1
	for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   

	    vector< const FstResult * > jacknivesToSend;
	    for(unsigned k=0;k<jackVec->size();k++){	       
		const FstResult * const* temsp= (jackVec->at(k)->getFstResult());
		//constFstResult
		//cout<<temsp[i][j]<<endl;
		//const FstResult * test1  = &temsp[i][j];
		//jacknivesToSend.push_back( &(  (jackVec->at(k)->getFstResult())[i][j]) );
		jacknivesToSend.push_back( &( temsp[i][j] ) );
	    }
	    toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<fstResults[i][j].printWithJacknife( &jacknivesToSend )<<endl;
	}
      }

      return toReturn.str();
    
}


void SumStatFst::computeStatSingle( const   AlleleRecords   * recordToUse,const bool allowUndefined){
    //currentRow = dataToUse->at(i);
    //cout<<"coord1 "<<dataToUse->at(d).coordinate<<" "<<d<<" "<<dataToUse->at(d) <<endl;
       
    if(!isResolvedDNA(recordToUse->ref) ){
	cerr<<"SumStateFst:  pairwiseFst() Problem for record #"<<"\t"<<recordToUse->chr<<" coordinate = "<<recordToUse->coordinate<<" reference = "<<recordToUse->ref<<" is not resolved"<<endl;
	exit(1);
    }


    //initialize
    if(!isResolvedDNA(recordToUse->alt)){ // if one of A,C,G,T
	return ; //next iteration, we need a valid alternative allele otherwise there is no impact on divergence calculation, this speeds it up
    }

    //taken from variant_file_output.cpp VCFtools
    int N_alleles = 2;
    vector<unsigned int> N_hom, N_het;
    vector<double> n(numberOfPopulations, 0.0);
    vector<vector<double> > p(numberOfPopulations, vector<double>(N_alleles,0.0));
    
    double nbar = 0.0;
    vector<double> pbar(N_alleles, 0.0);
    vector<double> hbar(N_alleles, 0.0);
    vector<double> ssqr(N_alleles, 0.0);
    double sum_nsqr = 0.0;
    double n_sum = 0.0;
    
    for (unsigned int i=2; i<numberOfPopulations; i++){
	//e->get_multiple_genotype_counts(indvs_in_pops[i], e->include_genotype, N_hom, N_het);
	

	if(!allowUndefined){
	    if(recordToUse->vectorAlleles->at(i).getTotalCount()<2){
		continue ;
	    }
	}
	double pRef             = double(recordToUse->vectorAlleles->at(i).getRefCount())/double(recordToUse->vectorAlleles->at(i).getTotalCount());
	double numDiploidInd = double(recordToUse->vectorAlleles->at(i).getTotalCount())/2.0;
	//j=0 ref
	N_hom[0] += pRef*pRef         * numDiploidInd;
	N_het[0] += pRef*(1-pRef)     * numDiploidInd;
		
	//j=0 alt
	N_hom[1] += (1-pRef)*(1-pRef) * numDiploidInd;
	N_het[1] += pRef*(1-pRef)     * numDiploidInd;

	for (unsigned int j=0; j<N_alleles; j++){
	    
	    cout<<"pop#"<<i<<" allele#"<<j<<"\thom"<<N_hom[j]<<"\thet"<<N_het[j]<<endl;
	
	
			    
	    n[i]   += N_hom[j] + 0.5*N_het[j];
	    p[i][j] = N_het[j] + 2*N_hom[j];
			    
	    nbar += n[i];
	    pbar[j] += p[i][j];
	    hbar[j] += N_het[j];
	}
	
	for (unsigned int j=0; j<N_alleles; j++)
	    p[i][j] /= (2.0*n[i]);	// diploid

	sum_nsqr += (n[i] * n[i]);
    }
	

#ifdef NO
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



    //the first is the ancestral, we should never reach that state given the check above
    if(sampledAllele[1] == 'N'){//the ancestral has an unresolved allele, skip
	//goto SKIPTONEXTITERATION;
	//continue;
	return ;
    }




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

	    computeFst(sampledAllele[1], //0 is root, 1 is ancestral
		       sampledAllele[i],
		       sampledAllele[j],
		       (cpgForPop[i] || cpgForPop[j]),
		       &(fstResults[i][j])  );
	       
	}
    }


#endif

    //SKIPTONEXTITERATION:
    return;
    
}



void SumStatFst::computeStat( const   vector < AlleleRecords >  * dataToUse,const vector<string> * popNames,const bool allowUndefined){
    numberOfPopulations=popNames->size()+1;//+1 for the human reference
    fstResults = new FstResult*[numberOfPopulations];
    for(unsigned int i=0;i<numberOfPopulations;i++)
	fstResults[i] = new FstResult[numberOfPopulations];

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
   // 	   cout<<i<<" "<<j<<" "<<fstResults[i][j]<<endl;
   ///++)
   //   delete(populationNames);

}


string SumStatFst::print() const {
    // cout<<"print() begin"<<endl;
    stringstream toReturn;
    for(unsigned i=2;i<numberOfPopulations;i++){
       for(unsigned j=2;j<numberOfPopulations;j++){	       
	   //skip when the population is identical
	   if(i==j)
	       continue;	   
	   toReturn<<populationNames->at(i)<<"-"<<populationNames->at(j)<<"\t"<<fstResults[i][j]<<endl;
       }
   }
    // cout<<"print() done"<<endl;
    return toReturn.str();
}

