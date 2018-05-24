/*
 * GlactoolsOperations
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlactoolsOperations.h"




//reading fasta index
void readFastaIndex(const string fastaIndex,
		    vector<chrinfo> & chrFound,
		    uint64_t & genomeLength){
    ifstream myFile;
    string line;

    myFile.open(fastaIndex.c_str(), ios::in);
    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    vector<string> fields = allTokens(line,'\t');
	    chrinfo toadd;

	    toadd.name         =fields[0];
	    toadd.startIndexChr=genomeLength+1;
	    toadd.length       =string2uint(fields[1]);
	    toadd.endIndexChr  =genomeLength+toadd.length;
	    chrFound.push_back(toadd);
	    genomeLength+=toadd.length;
	    // cout<<"gl "<<genomeLength<<endl;
	}
    }else{
	cerr<<"Cannot open fasta index  "<<fastaIndex<<endl;
	exit(1);
    }
    myFile.close();

}

//returns new defline
string initFiles(vector<GlacParser * > & vectorOfGP,
		 //bool & atLeastOneHasData,
		 vector<bool> & hasData,
		 vector<int> & popSizePerFile,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 //string & chr1,
		 uint16_t & chr1,
		 unsigned int & coordCurrent,
		 bool printOnlyFirstPop){

    string deflineToReturn = "#chr\tcoord\tREF,ALT\troot\tanc\t";


    //   atLeastOneHasData=false;
   hasData = vector<bool>(vectorOfGP.size(),false);


   for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	hasData[i] = vectorOfGP[i]->hasData()  ;
	if(!hasData[i]){
	    cerr<<"File #"<<(i+1)<<" does not have any data, exiting"<<endl;
	    exit(1);    
	}// else{
	//     //atLeastOneHasData=true;
	// }
    }


    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	unsigned int nonPop = vectorOfGP[i]->getPopulationsNames()->size();
	vector<string> pops;
	popSizePerFile.push_back(nonPop);

	for(unsigned int j=2;j<nonPop;j++){
	    pops.push_back( vectorOfGP[i]->getPopulationsNames()->at(j));
	}
	//cout<<vectorToString(pops,"\t");

	deflineToReturn+=vectorToString(pops,"\t");

	if(printOnlyFirstPop)
	    break;

	if( i!=(vectorOfGP.size() -1)){
	    if(!pops.empty())//if the current record was just root/anc
		deflineToReturn+="\t";
	    
	}       	
    }

    //deflineToReturn+="\n";



    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	vecAlleleRecords.push_back( vectorOfGP[i]->getData() );
	if(i==0){
	    chr1          = vecAlleleRecords[i]->chri;
	    coordCurrent  = vecAlleleRecords[i]->coordinate;
	}else{
	    if(chr1 != vecAlleleRecords[i]->chri ){
		cerr<<"initFiles() Chromosomes differ between "<<chr1<<" and "<< vecAlleleRecords[i]->chri <<endl;
		exit(1);    	
	    }
	    coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
	}	
    }
    return deflineToReturn;
}



bool sanityCheck(vector<GlacParser * > & vectorOfGP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 //string & chr1,
		 uint16_t & chr1,
		 unsigned int & coordCurrent,
		 //string & chrcheck ,
		 uint16_t & chrcheck,
		 char & refAllele,
		 bool force  ){
  
	    
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 

	if(hasData[i] && hasCoordinate[i]){
	    if(refAllele == '\0'){
		chrcheck   = vecAlleleRecords[i]->chri;
		refAllele  = vecAlleleRecords[i]->ref;
		if( chrcheck   != chr1){
		    cerr<<"sanityCheck()1 Chromosomes differ between "<<(* vecAlleleRecords[0])<<" and "<<(*vecAlleleRecords[i])<<endl;
		    exit(1);   
		}

	    }else{
		if( chrcheck   != vecAlleleRecords[i]->chri){
		    cerr<<"sanityCheck()2 Chromosomes differ between "<<(* vecAlleleRecords[0])<<" and "<<(*vecAlleleRecords[i])<<endl;
		    exit(1);   
		}

		if( refAllele  != vecAlleleRecords[i]->ref){
		    cerr<<"The reference allele differs between "<<(* vecAlleleRecords[0])<<" and "<<(*vecAlleleRecords[i])<<endl;
		    if(force)
			return false;
		    else
			exit(1);  
		}
	    }
	}
    }


    if(refAllele == '\0'){
	cerr<<"The reference allele could not be determined at coordinate "<<coordCurrent<<endl;	
	exit(1);  
    }

    return true;

}


bool printAllele(vector<GlacParser * > & vectorOfGP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<int> & popSizePerFile,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 uint16_t & chr1,
		 //string & chr1,
		 unsigned int & coordCurrent,
		 GlacWriter * gw,
		 bool isGL,
		 bool force){

    //sanity checks
    uint16_t chrcheck = 0;  //= vecAlleleRecords[0]->chr;
    char refAllele  = '\0'; // = vecAlleleRecords[0]->ref;
    bool isSane=sanityCheck(vectorOfGP,hasData,hasCoordinate,vecAlleleRecords,chr1,coordCurrent,chrcheck,refAllele,force);
    if(!isSane)
	return false;
    
    //ancestral info
    SingleAllele chimpAC;
    SingleAllele ancAC;
    SingleGL     chimpGL;
    SingleGL     ancGL;
    
    bool chimpSet=false;
    bool ancSet  =false;
	    
    //determining new alternative allele
    char newAlt = 'N';
    vector<SingleAllele> toPrint;
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){

	    if( !isResolvedDNA(newAlt)  && //is 'N'
		isResolvedDNA(vecAlleleRecords[i]->alt) ){ //not 'N'
		newAlt = vecAlleleRecords[i]->alt;
	    }
		
	    if( isResolvedDNA(newAlt)                   && //not 'N'
		isResolvedDNA(vecAlleleRecords[i]->alt) && //not 'N'
		vecAlleleRecords[i]->alt != newAlt){       //differ
		//tri-allelic
		//goto seekdata;
		return false;
	    }
	}
    }
	    
	  
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){
	    //chimp
	    //so ugly..
	    if(isGL){
		if(!vecAlleleRecords[i]->vectorGLs->at(0).alleleCountNull()){
		    if(chimpSet){			
			if(chimpGL != vecAlleleRecords[i]->vectorGLs->at(0)){
			    cerr<<"Disprency in chimp info between "<<(* vecAlleleRecords[i])<<" and "<<(chimpAC)<<endl;
			    if(force)
				return false;
			    else
				exit(1);  
			}
		    }else{
			chimpSet=true;
			chimpGL=vecAlleleRecords[i]->vectorGLs->at(0);
		    }
		}
	    }else{
		if(!vecAlleleRecords[i]->vectorAlleles->at(0).alleleCountNull()){
		    if(chimpSet){
			if(chimpAC != vecAlleleRecords[i]->vectorAlleles->at(0)){
			    cerr<<"Disprency in chimp info between "<<(* vecAlleleRecords[i])<<" and "<<(chimpAC)<<endl;
			    if(force)
				return false;
			    else
				exit(1);  
			}						
		    }else{
			chimpSet=true;
			chimpAC=vecAlleleRecords[i]->vectorAlleles->at(0);
		    }
		}
	    }

	    //anc
	    if(isGL){
		if(!vecAlleleRecords[i]->vectorGLs->at(1).alleleCountNull()){
		    if(ancSet){			
			if(ancGL != vecAlleleRecords[i]->vectorGLs->at(1)){
			    cerr<<"Disprency in ancestral info between "<<(* vecAlleleRecords[i])<<" and "<<(ancGL)<<endl;
			    if(force)
				return false;
			    else
				exit(1);  
			}
		    }else{
			ancSet=true;
			ancGL=vecAlleleRecords[i]->vectorGLs->at(1);
		    }
		}
	    }else{
		if(!vecAlleleRecords[i]->vectorAlleles->at(1).alleleCountNull()){
		    if(ancSet){
			if(ancAC != vecAlleleRecords[i]->vectorAlleles->at(1)){
			    cerr<<"Disprency in ancestral info between "<<(* vecAlleleRecords[i])<<" and "<<(ancAC)<<endl;
			    if(force)
				return false;
			    else
				exit(1);  
			}
			
		    }else{
			ancSet=true;
			ancAC=vecAlleleRecords[i]->vectorAlleles->at(1);
		    }
		}
	    }
	}
    }
    //cout<<endl;
	     

    // 	printnewrecord:
    //cout<<chr1<<"\t"<<coordCurrent<<"\t"<<refAllele<<","<<newAlt<<endl;
    //cout<<chr1<<"\t"<<coordCurrent<<"\t"<<refAllele<<","<<newAlt<<"\t"<<chimp<<"\t"<<anc<<"\t";
    AlleleRecords arw (isGL);


	
    
    arw.chri          = chr1;
    arw.coordinate    = coordCurrent;
    
    arw.ref           = refAllele;
    arw.alt           = newAlt;

    // cout<<"ref "<<arw.ref<<" >"<<refAllele<<"<"<<endl;
    // cout<<"alt "<<arw.alt<<" >"<<newAlt<<"<"<<endl;

    if(isGL){//ugly, is there a better way? inheritance?
	arw.vectorGLs = new vector<SingleGL>();
	arw.vectorGLs->push_back(chimpGL);
	arw.vectorGLs->push_back(ancGL);

    }else{
	arw.vectorAlleles = new vector<SingleAllele>();
	arw.vectorAlleles->push_back(chimpAC);
	arw.vectorAlleles->push_back(ancAC);
    }
    uint32_t sizePops=0;
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){
	    for(int k=2;k<popSizePerFile[i];k++){
		sizePops++;
		if(isGL)
		    arw.vectorGLs->push_back(vecAlleleRecords[i]->vectorGLs->at(k));
		    //toPrintGL.push_back(vecAlleleRecords[i]->vectorGLs->at(k));
		else
		    arw.vectorAlleles->push_back(vecAlleleRecords[i]->vectorAlleles->at(k));
		    //toPrintAC.push_back(vecAlleleRecords[i]->vectorAlleles->at(k));
	    }
		    
	}else{

	    for(int k=2;k<popSizePerFile[i];k++){
		sizePops++;
		if(isGL){
		    SingleGL t;
		    arw.vectorGLs->push_back(t);
		}else{
		    SingleAllele t;
		    ///toPrintAC.push_back(t);
		    arw.vectorAlleles->push_back(t);
		}
	    }

	}
    }
    arw.sizePops      = sizePops;

    //cout<<vectorToString(toPrint,"\t")<<endl;
    if(!gw->writeAlleleRecord(&arw)){
	cerr<<"GlactoosOperation: printAllele() error record "<<arw<<endl;
	exit(1);
    }

	//co
    return true;
}




//! Method to read a sorted bed file (copy pasted from vcfcompute.cpp (did not integrate well with utils.h))
/*!
 *
 * This method checks for the records being ordered

  \param filetoread : String with the full path to the file to read
  \return           : Return(head) a pointer to a map where the key is the chromosome name and the value a vector of genomic ranges
  \sa  readBEDSortedfile()
*/

map< string, vector<GenomicRange> * > * readBEDSortedfile(string filetoread){
    //vector<GenomicRange> toReturn;
    map< string, vector<GenomicRange> * > * toReturn= new map< string, vector<GenomicRange> *>();
    igzstream myFile;
    myFile.open(filetoread.c_str(), ios::in);
    string line;
    unsigned int     lastEndCoord = 0;
    string           lastChrname  = "###";

    if (myFile.good()){
	while ( getline (myFile,line)){	    
	    //cout<<line<<endl;
	    if(line.empty()){
		continue;
	    }
	    string       chrName;
	    unsigned int startCoord;
	    unsigned int endCoord;
	    vector<string> temp=allTokens(line,'\t');
	    if(temp.size() != 3){
		cerr << "Error in readBEDSortedfile(): following line does not have 3 fields"<< line<<endl;
		exit(1);		
	    }
	    chrName     = destringify<string>(temp[0]);
	    startCoord  = destringify<unsigned int>(temp[1])+1; //the left coordinate is zero based
	    endCoord    = destringify<unsigned int>(temp[2]);

	    if(lastChrname != chrName){//new chr found
		//check for previously found		
		lastChrname  = chrName;
		lastEndCoord = endCoord;
		if(toReturn->find(chrName) == toReturn->end() ){//not previously found
		    //cout<<"new chr "<<chrName<<endl;
		    (*toReturn)[chrName]=new vector<GenomicRange>();
		}else{
		    cerr << "Cannot have multiple records on the same chromosome at different parts of the file "<< line<<endl;
		    exit(1);
		}
	    }else{//stay on same chr
		if(startCoord < lastEndCoord ){
		    cerr << "Problem with line =  "<<line<<" the start of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
		    exit(1);
		}

		if(endCoord   < lastEndCoord ){
		    cerr << "Problem with line =  "<<line<<" the end of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
		    exit(1);		    
		}
	    }

	    GenomicRange toadd (chrName,startCoord,endCoord);	    
	    (*toReturn)[chrName]->push_back(toadd);
	    //cout<<(*toReturn)[chrName]->size()<<endl;
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filetoread<<endl;
	exit(1);
    }

    return toReturn;
}


