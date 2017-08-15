/*
 * VCFreader
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "VCFreader.h"
// #define DEBUG


VCFreader::VCFreader(string file,string indexForFile,string chrName,int start,int end,int indelsAhead){
    readAhead=indelsAhead;
    rt =new ReadTabix (file,indexForFile,chrName,start,end);

    needToPopulateQueue=true;
    fullQueue          =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;

    svcfToReturn=0;

    repoCalledHasData=false;
	    
    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    tabixMode = true;
    textMode  = false;
}




VCFreader::VCFreader(string file,int indelsAhead){
    readAhead=indelsAhead;
    numberOfTimesHasDataWasCalled=0;
    svcfToReturn=0;

    vcfFile.open(file.c_str(), ios::in);    // open the streams
    if (vcfFile.good()) {
	//fine
    }else{
	cerr<<"Unable to open the file "<<file<<endl;
	exit(1);

    }

    needToPopulateQueue =true;
    fullQueue           =false;
    endQueue            =false;
    repoCalledHasData=false;


    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    tabixMode = false;
    textMode  = true;
}




VCFreader::~VCFreader(){

    if(tabixMode){
	delete rt; //calling the destructor
    }

    if(textMode)
	vcfFile.close();

    // if(!queueOfVCFs.empty()){
    // 	cerr<<"The queue still contains elements " <<endl;
    // 	exit(1);
    // }

    while(!queueOfVCFs.empty()){
#ifdef DEBUG			
	cout<<"deleteQueue "<<svcfToReturn<<endl;
#endif
	delete( queueOfVCFs.front() );
	queueOfVCFs.pop_front();
    }

    queueOfVCFs.clear();
#ifdef DEBUG		
    cout<<"deleteDestr "<<svcfToReturn<<endl;
#endif
    delete svcfToReturn;
}

void VCFreader::repositionIterator(string chrName,int start,int end){


    if(!tabixMode){
	cerr<<"The subroutine repositionIterator can only be called on objects constructed using tabix " <<endl;
	exit(1);	
    }


    //re-initializing variables
    needToPopulateQueue =true;
    fullQueue           =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;
#ifdef DEBUG			
    cout<<"deleteRepo "<<svcfToReturn<<endl;
#endif
    delete svcfToReturn;
    svcfToReturn=0;

    
    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    while(!queueOfVCFs.empty()){
#ifdef DEBUG			
    cout<<"deleteQue "<<svcfToReturn<<endl;
#endif
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
	SimpleVCF * current=queueOfVCFs.front();
	//SimpleVCF * current=getData();
	//cout<<current->getPosition()<<endl;
	//when the current position is found, we set the flag repoCalledHasData
	//such that hasData will return true and the first element of the queue will be returned
	if(int(current->getPosition()) >= start){
	    repoCalledHasData=true;
	    break;
	}
	current=getData();
    }
    

}

bool VCFreader::getNextLine(){
    if(tabixMode){
	return rt->readLine(currentline);
    }

    if(textMode){
	while(1){
	    bool flag=(bool)getline(vcfFile,currentline);
	    if(!flag)
		return false;	   
	    if(currentline.length() > 0 && currentline[0] != '#')
		return true;
	}
    }
    cerr<<"Invalid state in VCFreader::getNextLine()"<<endl;
    exit(1);
    return false;
}


bool VCFreader::hasData(){

    if(repoCalledHasData){
	repoCalledHasData=false;
	return true;
    }

    numberOfTimesHasDataWasCalled++;

    
    //if first call and queue empty, populate
    if(needToPopulateQueue){
	bool loop=true;
	int indexQueue=0;
	while(loop){
	    if(getNextLine()){
#ifdef DEBUG		
		cout<<"currentline "<<currentline<<endl;
#endif
		SimpleVCF * svcf = new SimpleVCF(currentline);
#ifdef DEBUG		
		cout<<"new1 "<<svcf<<endl;
#endif
		if(queueOfVCFs.size() != 0 ){
		    flagCpG( queueOfVCFs.back(),svcf);
		}
		// cout<<"Adding "<<*svcf<<endl;
		queueOfVCFs.push_back(svcf);

		if(svcf->containsIndel()){
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
	    SimpleVCF * svcf = new  SimpleVCF(currentline);
#ifdef DEBUG		
	    cout<<"new2 "<<svcf<<endl;
	    cout<<"size "<<queueOfVCFs.size()<<endl;
#endif

	    if(queueOfVCFs.size() != 0 ){
		flagCpG( queueOfVCFs.back(),svcf);
	    }	    
	    queueOfVCFs.push_back(svcf);
	    
	    if(svcf->containsIndel()){	
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


inline void VCFreader::flagCpG(SimpleVCF * previous,SimpleVCF * current){ //pass by address
    //cout<<"flagCpG "<<previous->getPosition()<<"\t"<<previous->hasAtLeastOneC()<<"\t"<<current->getPosition()<<"\t"<<current->hasAtLeastOneG()<<endl;
    if( ( (previous->getPosition()+1) == current->getPosition())  &&   //one position behind
	(  previous->getChr()         == current->getChr()    )   &&   //on same chr
	(previous->hasAtLeastOneC()   && current->hasAtLeastOneG()) ){  //previous has at least one C, current has at least one G
	previous->setCpg(true);
	current->setCpg(true);
    }
}



SimpleVCF * VCFreader::getData(){
// auto_ptr<AlleleInfo> VCFreader::getData(){

    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;

    //delete the previous data
#ifdef DEBUG			
    cout<<"delete2 "<<svcfToReturn<<endl;
#endif

    delete svcfToReturn;
    //}

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
	    svcfToReturn->setCloseIndel(true);
	    indexInQueueOfIndels = -1;
	}else{
	    if(indexInQueueOfIndels != 0){//elements in the list that are indels, need to check them
		int indexInList=0;	    
		list<SimpleVCF *>::iterator it;
		for ( it=queueOfVCFs.begin() ; it != queueOfVCFs.end(); it++ ){
		    if(indexInList<=indexInQueueOfIndels){
			if( (*it)->containsIndel()){

			    if( ( int((*it)->getPosition() - svcfToReturn->getPosition() ) <= readAhead) &&
				((*it)->getChr() == svcfToReturn->getChr()) ){
				svcfToReturn->setCloseIndel(true);			    
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
	 if( (int(svcfToReturn->getPosition() - indexOfLastIndel) <= readAhead) &&
	     (chrOfLastIndel == svcfToReturn->getChr()) ){
	     svcfToReturn->setCloseIndel(true);
	 }
     }

     

    if(svcfToReturn->containsIndel()){
	previouslyFoundIndel=true;
	indexOfLastIndel    =svcfToReturn->getPosition();
	chrOfLastIndel      =svcfToReturn->getChr();
    }

    return svcfToReturn;
}
