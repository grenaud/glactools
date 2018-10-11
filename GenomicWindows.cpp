/*
 * GenomicWindows
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GenomicWindows.h"


GenomicWindows::GenomicWindows(string sqLines,bool allowSexChr){
    //    srand48 ( time(NULL) );
    // time_t t; 
    // (void) time(&t);
    // cout<<"t "<<t<<endl;

    timeval time;
    gettimeofday(&time, NULL);
    srand48(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );
    genomeLength=0;
    this->allowSexChr=allowSexChr;

    vector<string> sq = allTokens(sqLines,'\n');
    for(unsigned int i=0;i<sq.size();i++){
    	vector<string> sql = allTokens(sq[i],'\t');
    	//cout<<sql[2]<<endl;	
	if(sql.size()!=3){
	    cerr<<"ERROR line from SQ header does not have 3 fields"<<sq[i]<<endl;
	}
	chrinfo toadd;
	    
	toadd.name         = sql[1].substr(3);
	toadd.startIndexChr=genomeLength+1;
	toadd.length       = string2uint(sql[2].substr(3));//string2uint(fields[1]);
	toadd.endIndexChr  =genomeLength+toadd.length;
	chrFound.push_back(toadd);
	genomeLength+=toadd.length;	    
    }
	


    //BEGIN Reading the fasta index 
    // ifstream myFile;
    // string line;
    // myFile.open(fastaIndex.c_str(), ios::in);
    // if (myFile.is_open()){
    // 	while ( getline (myFile,line)){
    // 	    vector<string> fields = allTokens(line,'\t');
    // 	    chrinfo toadd;

    // 	    toadd.name         =fields[0];
    // 	    toadd.startIndexChr=genomeLength+1;
    // 	    toadd.length       =string2uint(fields[1]);
    // 	    //cout<< toadd.name <<"\t"<<toadd.length<<endl;
    // 	    toadd.endIndexChr  =genomeLength+toadd.length;
    // 	    chrFound.push_back(toadd);
    // 	    genomeLength+=toadd.length;
    // 	}
    // }else{
    // 	cerr<<"Cannot open fasta index  "<<fastaIndex<<endl;
    // 	exit(1);
    // }
    // myFile.close();
    //END Reading the fasta index

}

GenomicWindows::~GenomicWindows(){

}






vector<GenomicRange> GenomicWindows::getGenomeWide(){

    vector<GenomicRange> toReturn;

    for(unsigned int i=0;i<chrFound.size();i++){ 
	vector<GenomicRange> toAdd=getChr(chrFound[i].name);
	if(!toAdd.empty()){
	    for(unsigned int k=0;k<toAdd.size();k++){
		toReturn.push_back(toAdd[k]);
	    }
	}
    }



    return toReturn;
}

vector<GenomicRange> GenomicWindows::getChr(string chrName){
    vector<GenomicRange> toReturn;

    if(!allowSexChr){
	if( (chrName.find("X") != string::npos) ||
	    (chrName.find("Y") != string::npos) ){
	    return toReturn;//empty vector	    
	}
    }

    bool foundChr=false;
    int indexChr=-1;
    for(unsigned int i=0;i<chrFound.size();i++){
	if(chrFound[i].name == chrName){
	    foundChr=true;
	    indexChr=i;
	}
    }


    if(!foundChr){
	cerr<<"Error: Cannot find chromosome "<<chrName<<" exiting"<<endl;
	exit(1);
    }


    GenomicRange toAdd;
    toAdd.setChrName( chrFound[indexChr].name);
    toAdd.setStartCoord( 1);
    toAdd.setEndCoord( chrFound[indexChr].length);
    toReturn.push_back(toAdd);

    return toReturn;

}






vector<GenomicRange> GenomicWindows::getGenomicWindowsChr(string chrName,int windowSize,int overlap){
    if(overlap >= windowSize){
	cerr<<"Cannot have an overlap "<<overlap<<" greater than the "<<windowSize<<endl;
	exit(1);
    }

    vector<GenomicRange> toReturn;

    if(!allowSexChr){
	if( (chrName.find("X") != string::npos) ||
	    (chrName.find("Y") != string::npos) ){
	    return toReturn;//empty vector	    
	}
    }

    bool foundChr=false;
    int indexChr=-1;
    for(unsigned int i=0;i<chrFound.size();i++){
	if(chrFound[i].name == chrName){
	    foundChr=true;
	    indexChr=i;
	}
    }

    if(!foundChr){
	cerr<<"Error: Cannot find chromosome "<<chrName<<" exiting"<<endl;
	exit(1);
    }

    int coord =1;
    while(//chrFound[indexChr].startIndexChr <= coord &&	  
	  coord                     <= (int(chrFound[indexChr].length)-windowSize)){
	GenomicRange toAdd;
	toAdd.setChrName( chrFound[indexChr].name);
	toAdd.setStartCoord( coord);
	toAdd.setEndCoord(   coord+windowSize-1);
	toReturn.push_back(toAdd);
	coord+=(windowSize-overlap);
    }

    return toReturn;
}

vector<GenomicRange> GenomicWindows::getGenomicWindows(int windowSize,int overlap){

    if(overlap >= windowSize){
    	cerr<<"Cannot have an overlap "<<overlap<<" greater than the "<<windowSize<<endl;
    	exit(1);
    }

    vector<GenomicRange> toReturn;


    for(unsigned int i=0;i<chrFound.size();i++){ 
	//toReturn.push_back(getGenomicWindowsChr(chrFound[i].name,windowSize,overlap));
	vector<GenomicRange> toAdd=getGenomicWindowsChr(chrFound[i].name,windowSize,overlap);
	//	cout<<chrFound[i].name<<"\t"<<toAdd.size()<<endl;
	if(!toAdd.empty()){
	    for(unsigned int k=0;k<toAdd.size();k++){
		//cout<<toAdd[k]<<endl;
		toReturn.push_back(toAdd[k]);
	    }
	}
    }



    return toReturn;
}
