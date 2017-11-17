/*
 * MultiVCFreader
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "MultiVCFreader.h"
//#define DEBUG

//IMPLEMENT HEADER
MultiVCFreader::MultiVCFreader(string file,string indexForFile,string chrName,int start,int end,int indelsAhead){
    readAhead=indelsAhead;
    rt = new ReadTabix (file,indexForFile,chrName,start,end);
    strTbx = rt->getKstringPtr();
    str = strTbx;
    memset(&aux, 0, sizeof(ks_tokaux_t));
    
    //    cerr<<"first const"<<endl;
    //reading header
    istringstream f (rt->getHeader());
    string line;    
    numPop=0;
    // V not exactly optimized but this part is not time critical V
    while (getline(f, line)) {
        //std::cout << line << std::endl;
	if(strBeginsWith(line,"#CHROM") 
	   // ||
	   // strBeginsWith(line,"#chr")
	){
	    
	    vector<string> tok = allTokens(line,'\t');
	    if(tok.size() < 10 ){
		cerr<<"The header line"<<line<<" does not contain enough fields"<<endl;
		exit(1);		
	    }

	    for(unsigned int i=9;i<tok.size();i++){
		//cerr<<tok[i]<<endl;
		numPop++;
		populationNames.push_back(tok[i]);
	    }

	}
    }


    if( numPop == 0 ){
	cerr<<"No populations have been found for file "<<file<<endl;
	exit(1);
    }

    needToPopulateQueue           = true;
    fullQueue                     = false;
    endQueue                      = false;
    numberOfTimesHasDataWasCalled = 0;

    svcfToReturn                  = 0;

    repoCalledHasData             = false;
	    
    indexInQueueOfIndels          = -1;
    indexOfLastIndel              = 0;
    previouslyFoundIndel          = false;

    tabixMode                     = true;
    textMode                      = false;
}




MultiVCFreader::MultiVCFreader(string file,int indelsAhead){
    readAhead=indelsAhead;
    numberOfTimesHasDataWasCalled=0;
    svcfToReturn=0;

    //vcfFile.open(file.c_str(), ios::in);    // open the streams

    if( (fpVCF = bgzf_open(file.c_str(), "r")) == 0){ 
    	cerr<<"Unable to open the file "<<file<<endl;
    	exit(1);
    }

    // if (vcfFile.good()) {
    // 	//fine
    // }else{
    // 	cerr<<"Unable to open the file "<<file<<endl;
    // 	exit(1);
    // }

    bool firstLine=true;
    bool haveCaptureCHROM=false;
    numPop=0;


    while(1){
	//bool flag=(bool)getline(vcfFile,currentline);
	int ret;
	strText= new kstring_t();
	strText->s = 0; strText->l = strText->m = 0;
	str = strText;
        memset(&aux, 0, sizeof(ks_tokaux_t));

	bool flag= ((ret = bgzf_getline(fpVCF, '\n', strText)) >= 0); 
	string currentline (strText->s); //copy but this is just the header

	if(!flag){
	    cerr<<"ERROR file : "+file+" is probably empty"<<endl;
	    exit(1);
	}

	if(firstLine){

	    if(currentline.length() > 0 && currentline[0] != '#'){
	    
		cerr<<"ERROR first line in "<<file<<"does not begin with #"<<endl;
		exit(1);
	    }

	    firstLine = false;
	}

	if(!firstLine){

	    
	    if(currentline.length() > 0){
		if(currentline[0] == '#'){


		    if(strBeginsWith(currentline,"#CHROM") 
		       // ||
		       // strBeginsWith(currentline,"#chr") 
		    ){
			haveCaptureCHROM=true;
			vector<string> tok = allTokens(currentline,'\t');
			if(tok.size() < 10 ){
			    cerr<<"The header line"<<currentline<<" does not contain enough fields for file "<<file<<endl;
			    exit(1);		
			}
			
			for(unsigned int i=9;i<tok.size();i++){
			    //cerr<<tok[i]<<endl;
			    numPop++;
			    populationNames.push_back(tok[i]);
			}
			break;
		    }

		    // cerr<<"ERROR first line in "<<file<<"does not begin with #"<<endl;
		    // return 1;
		}else{
		    break;
		}
	    }

	}

    }//end while(1)
    // vcfFile.close();
    // //vcfFile.seekg(0, std::ios::beg);

    // vcfFile.open(file.c_str(), ios::in);    // open the streams
    // if (vcfFile.good()) {
    // 	//fine
    // }else{
    // 	cerr<<"Unable to open the file for second pass "<<file<<endl;
    // 	exit(1);
    // }

    if( numPop == 0 ){
	cerr<<"No populations have been found for file "<<file<<endl;
	exit(1);
    }

    if(!haveCaptureCHROM){
	cerr<<"The header with #CHROM has not been found in file:"<<file<<endl;
	exit(1);		
    }
    

    needToPopulateQueue = true;
    fullQueue           = false;
    endQueue            = false;				   
    repoCalledHasData   = false;


    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    tabixMode = false;
    textMode  = true;
}




MultiVCFreader::~MultiVCFreader(){

    if(tabixMode){
	delete rt; //calling the destructor
    }

    if(textMode)
	delete strText;
    //vcfFile.close();
    bgzf_close(fpVCF);
    


    // if(!queueOfVCFs.empty()){
    // 	cerr<<"The queue still contains elements " <<endl;
    // 	exit(1);
    // }

    while(!queueOfVCFs.empty()){
	delete( queueOfVCFs.front() );
	queueOfVCFs.pop_front();
    }

    queueOfVCFs.clear();

    //cout<<"deleteDest "<<svcfToReturn<<endl;
    //cout<<"delete3 "<<svcfToReturn->size()<<endl;
    if(svcfToReturn != 0){
	for(unsigned int i=0;i<svcfToReturn->size();i++){
	    delete svcfToReturn->at(i);
	}
	delete svcfToReturn;
    }
}

void MultiVCFreader::repositionIterator(string chrName,int start,int end){


    if(!tabixMode){
	cerr<<"The subroutine repositionIterator can only be called on objects constructed using tabix " <<endl;
	exit(1);	
    }


    //re-initializing variables
    needToPopulateQueue =true;
    fullQueue           =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;

    //cout<<"delete1 "<<svcfToReturn->size()<<endl;
    if(svcfToReturn != 0){
    for(unsigned int i=0;i<svcfToReturn->size();i++){
	delete svcfToReturn->at(i);
    }
    delete svcfToReturn;
    }
    svcfToReturn=0;

    
    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    while(!queueOfVCFs.empty()){
	delete( queueOfVCFs.front() );
	queueOfVCFs.pop_front();
    }

    queueOfVCFs.clear();

    //if(bp
    //we need to jump a bit behind to detect CpG and indels
    int coordinatePrior=max(readAhead,1);
    rt->repositionIterator(chrName,start-coordinatePrior,end);
    

    //need to reposition the queue to the desired coord
    while(hasData()){
	vector<SimpleVCF *> * current=queueOfVCFs.front();

	//SimpleVCF * current=getData();
	//cout<<current->getPosition()<<endl;
	//when the current position is found, we set the flag repoCalledHasData
	//such that hasData will return true and the first element of the queue will be returned
	if(int(current->at(0)->getPosition()) >= start){
	    repoCalledHasData=true;
	    break;
	}
	current=getMultipleData();
    }
    

}

vector<string> MultiVCFreader::getPopulationNames() const {
    return populationNames;
}

bool MultiVCFreader::getNextLine(){
    if(tabixMode){
	//return rt->readLine(currentline);
	return rt->readLineKS();//string is found  in strTbx
    }

    if(textMode){
	while(1){
	    //bool flag=(bool)getline(vcfFile,currentline);
	    int ret;
	    bool flag= ((ret = bgzf_getline(fpVCF, '\n', strText)) >= 0); 

	    if(!flag)
		return false;	   
	    //if(currentline.length() > 0 && currentline[0] != '#')
	    if(strText->l > 0 && strText->s[0] != '#')
		return true;
	}
    }

    cerr<<"Invalid state in MultiVCFreader::getNextLine()"<<endl;
    exit(1);

    return false;
}


bool MultiVCFreader::hasData(){
    //cout<<"hasData"<<endl;
    if(repoCalledHasData){
	repoCalledHasData=false;
	return true;
    }

    numberOfTimesHasDataWasCalled++;

    
    //if first call and queue empty, populate
    if(needToPopulateQueue){
#ifdef DEBUG		
	cerr<<"MultiVCFreader hasData()"<<endl;
#endif
	bool loop=true;
	int indexQueue=0;
	while(loop){
	    if(getNextLine()){
#ifdef DEBUG		
		cerr<<"currentline "<<currentline<<endl;
#endif
		vector<SimpleVCF *> *  svcfvec = new vector<SimpleVCF *>();
		//vector<string> fieldTab = allTokens(currentline,'\t');
		char *p, *q;
		int i;
		CoreVCF * corevcf =  new CoreVCF();
		int k=0;
		for (p = kstrtok(str->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
		    q = (char*)aux.p;
		    *q = 0;

		    if(i<9){//core VCF fields
			if(i == 0){ corevcf->setName(  p); continue; }
			if(i == 1){ corevcf->setPos(   p); continue; }
			if(i == 2){ corevcf->setID(    p); continue; }
			if(i == 3){ corevcf->setREF(   p); continue; }
			if(i == 4){ corevcf->setALT(   p); continue; }
			if(i == 5){ corevcf->setQUAL(  p); continue; }
			if(i == 6){ corevcf->setFILTER(p); continue; }
			if(i == 7){ corevcf->setINFO(  p); continue; }
			if(i == 8){ corevcf->setFORMAT(p); continue; }
		    }else{
			
			SimpleVCF * svcf = new  SimpleVCF (p,corevcf,k==0);
			svcfvec->push_back(svcf);
			k++;
		    }

		}

		if( (9+numPop) != i){
		    cerr<<"Found "<<i<<" fields but expected "<<(9+numPop)<<endl;
		    exit(1);
		}
		// 		for(int k=0;k<numPop;k++){
     
// 		    SimpleVCF * svcf = new  SimpleVCF (fieldTab,corevcf,k==0);
		    
// #ifdef DEBUG		
// 		    cerr<<"field#"<<k<<"/"<<numPop<<" = ->"<<fieldTab[k+9]<<"<- pos="<<svcf->getPosition()<<endl;
// #endif

// 		    svcfvec->push_back(svcf);
// 		}
		
		//SimpleVCF * svcfvec = new SimpleVCF(currentline);
#ifdef DEBUG		
		//cerr<<"new1 "<<svcf<<endl;
		cerr<<"done "<<endl;
#endif
		
		if(queueOfVCFs.size() != 0 ){
		    //for(
		    for(int k=0;k<numPop;k++)
			flagCpG( queueOfVCFs.back()->at(k) , svcfvec->at(k) );
		}
#ifdef DEBUG
	       cerr<<"Adding "<<endl;
#endif
		queueOfVCFs.push_back(svcfvec);

		if(svcfvec->at(0)->containsIndel()){
		    if(indexInQueueOfIndels == -1)
			indexInQueueOfIndels=indexQueue;		   
		}
		indexQueue++;
		if(int(queueOfVCFs.size()) == (readAhead+1)){
		    loop=false;
		}
	    }else{
		loop=false;
	    }
	}
	

	if(int(queueOfVCFs.size()) == (readAhead+1)){ //+1 for CPGs
	    fullQueue=true;
	}else{
	    endQueue=true;
	}

	needToPopulateQueue=false;
    }

    //if subsequent call, and queue full
    if(fullQueue){
	if(getNextLine()){
	    // SimpleVCF * svcf = new  SimpleVCF(currentline);

	    vector<SimpleVCF *> *  svcfvec = new vector<SimpleVCF *>();

	    char *p, *q;
	    int i;
	    CoreVCF * corevcf =  new CoreVCF();
	    int k=0;
	    for (p = kstrtok(str->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
		q = (char*)aux.p;
		*q = 0;

		//cerr<<"tokenizer token#1"<<i<<"="<<p<<"#"<<endl;
		if(i<9){//core VCF fields
		    if(i == 0){ corevcf->setName(  p); continue; }
		    if(i == 1){ corevcf->setPos(   p); continue; }
		    if(i == 2){ corevcf->setID(    p); continue; }
		    if(i == 3){ corevcf->setREF(   p); continue; }
		    if(i == 4){ corevcf->setALT(   p); continue; }
		    if(i == 5){ corevcf->setQUAL(  p); continue; }
		    if(i == 6){ corevcf->setFILTER(p); continue; }
		    if(i == 7){ corevcf->setINFO(  p); continue; }
		    if(i == 8){ corevcf->setFORMAT(p); continue; }
		}else{
			
		    SimpleVCF * svcf = new  SimpleVCF (p,corevcf,k==0);
		    svcfvec->push_back(svcf);
		    k++;
		}

	    }

	    if( (9+numPop) != i){
		cerr<<"Found "<<i<<" fields but expected "<<(9+numPop)<<endl;
		exit(1);
	    }

	    // vector<string> fieldTab = allTokens(currentline,'\t');
	    // CoreVCF * corevcf =  new CoreVCF(fieldTab);
	    // for(int k=0;k<numPop;k++){
	    // 	SimpleVCF * svcf = new SimpleVCF (fieldTab,corevcf,k==0);
	    // 	svcfvec->push_back(svcf);
	    // }
		
#ifdef DEBUG		
	    //cout<<"new2 "<<*svcf<<endl;
#endif
	    // cout<<"size "<<queueOfVCFs.size()<<endl;
	    if(queueOfVCFs.size() != 0 ){
		for(int k=0;k<numPop;k++)
		    flagCpG( queueOfVCFs.back()->at(k),svcfvec->at(k));
	    }
	    
	    queueOfVCFs.push_back(svcfvec);
	    
	    if(svcfvec->at(0)->containsIndel()){	
		if(indexInQueueOfIndels == -1)
		  indexInQueueOfIndels=queueOfVCFs.size()-1;
	    }
	    
	}else{
	    fullQueue=false;
	    endQueue=true;
	}

    }

    //if final calls and queue not max size
    if(endQueue){
	//nothing to do
    }

    
    bool stillHasData=!(queueOfVCFs.empty());

    // if(!stillHasData){ //getData() should not get called in this case, hence no deallocation
    // 	cout<<"delete1  "<<svcfToReturn<<endl;
    // 	delete svcfToReturn;
    // }

    return stillHasData;
}


inline void MultiVCFreader::flagCpG(SimpleVCF * previous,SimpleVCF * current){ //pass by address

    // cerr<<"flagCPG"<<endl;
    // cerr<<(previous->getPosition()+1)<<endl;
    // cerr<<current->getPosition()<<endl;
    // cerr<<previous->getChr()<<endl;
    // cerr<<current->getChr()<<endl;
    // cerr<<previous->hasAtLeastOneC()<<endl;
    // cerr<<current->hasAtLeastOneG()<<endl;
    
    if( ( (previous->getPosition()+1) == current->getPosition())  &&   //one position behind
	(  previous->getChr()         == current->getChr()    )   &&   //on same chr
	(previous->hasAtLeastOneC()   && current->hasAtLeastOneG()) ){  //previous has at least one C, current has at least one G

	
	// exit(1);
	previous->setCpg(true);
	current->setCpg(true);
    }
}


SimpleVCF  * MultiVCFreader::getData(){
    cerr<<"Cannot call getData for a MultiVCFreader "<<endl;
    exit(1);
}

vector<SimpleVCF *> * MultiVCFreader::getMultipleData(){
// auto_ptr<AlleleInfo> MultiVCFreader::getData(){
    //cout<<"getMultipleData"<<endl;
    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;

    //delete the previous data

    //cout<<"delete2"<<endl;
    if(svcfToReturn != 0 ){
	//cout<<"delete2 "<<svcfToReturn->size()<<endl;
	for(unsigned int i=0;i<svcfToReturn->size();i++){
	    //cout<<"delete2 "<<i<<endl;
	    delete svcfToReturn->at(i);
	}
	delete svcfToReturn;
    }

    //SimpleVCF svcf =queueOfVCFs.front();

    //svcfToReturn =new SimpleVCF(queueOfVCFs.front());
    svcfToReturn = queueOfVCFs.front();

    // cout<<&(queueOfVCFs.front())<<endl;
    // cout<<svcfToReturn<<endl;
    //delete svcfToReturn;
    queueOfVCFs.pop_front();


    // cout<<endl<<"indexInQueueOfIndels "<<indexInQueueOfIndels<<endl;
    //look ahead   
    if(readAhead != 0){//we need to set the close to indel flag
	if(indexInQueueOfIndels == 0){//last element that is to be returned is an indel
	    svcfToReturn->at(0)->setCloseIndel(true);
	    indexInQueueOfIndels = -1;
	}else{
	    if(indexInQueueOfIndels != 0){//elements in the list that are indels, need to check them
		int indexInList=0;	    
		list< vector<SimpleVCF * > * >::iterator it;

		for ( it=queueOfVCFs.begin() ; it != queueOfVCFs.end(); it++ ){
		    if(indexInList <= indexInQueueOfIndels){
			//it->at(0)->getCorevcf()

			if( (*it)->at(0)->containsIndel()){

			    if( (  int((*it)->at(0)->getPosition() - svcfToReturn->at(0)->getPosition())  <= readAhead) &&
				((*it)->at(0)->getChr() == svcfToReturn->at(0)->getChr()) ){
				svcfToReturn->at(0)->setCloseIndel(true);			    
			    }

			}
		    }else{
			break;
		    }
		
		    indexInList++;
		}

	    }
	    indexInQueueOfIndels=max(indexInQueueOfIndels-1,-1);	
	}
    }

    //look behind for indels
    if(previouslyFoundIndel){
	if( ( int(svcfToReturn->at(0)->getPosition() - indexOfLastIndel) <= readAhead) &&
	    (chrOfLastIndel == svcfToReturn->at(0)->getChr()) ){
	    svcfToReturn->at(0)->setCloseIndel(true);
	}
     }

     

    if(svcfToReturn->at(0)->containsIndel()){
	previouslyFoundIndel=true;
	indexOfLastIndel    =svcfToReturn->at(0)->getPosition();
	chrOfLastIndel      =svcfToReturn->at(0)->getChr();
    }

    return svcfToReturn;
}
