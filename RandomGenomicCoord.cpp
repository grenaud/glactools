/*
 * RandomGenomicCoord
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "RandomGenomicCoord.h"



uint64_t   RandomGenomicCoord::randGenomicCoord_(){     
     while(true){
	 uint64_t toReturn = 
	     (((uint64_t) rand() <<  0) & 0x000000000000FFFFull) | 
	     (((uint64_t) rand() << 16) & 0x00000000FFFF0000ull) | 
	     (((uint64_t) rand() << 32) & 0x0000FFFF00000000ull) |
	     (((uint64_t) rand() << 48) & 0xFFFF000000000000ull);
	 //cout<<toReturn<<endl;
	 return toReturn%genomeLength;//should work if genome size is less than 2^32
	 //	 if(toReturn<genomeLength)
	 //	     return toReturn;
     }
 }

RandomGenomicCoord::RandomGenomicCoord(string linesSQ,bool allowSexChr){
    //    srand48 ( time(NULL) );
    // time_t t; 
    // (void) time(&t);
    // cout<<"t "<<t<<endl;

    timeval time;
    gettimeofday(&time, NULL);
    srand48(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );
    genomeLength=0;
    this->allowSexChr=allowSexChr;
    //BEGIN Reading the fasta index


    vector<string> sq = allTokens(linesSQ,'\n');
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
	

}

RandomGenomicCoord::~RandomGenomicCoord(){

}

GenomicRange RandomGenomicCoord::getRandomGenomic(int bpToExtract){
    bool found=false;
    GenomicRange toReturn;


    while(!found){
	uint64_t coord =randGenomicCoord_();

	for(unsigned int i=0;i<chrFound.size();i++){
	    
	    if(chrFound[i].startIndexChr <= coord &&
	       coord                     <= (chrFound[i].endIndexChr-bpToExtract)){
		found=true;


		toReturn.setChrName(       chrFound[i].name);
		toReturn.setStartCoord( coord-chrFound[i].startIndexChr);
		toReturn.setEndCoord(   coord-chrFound[i].startIndexChr+bpToExtract);
		
		if(!allowSexChr){
		    if( (toReturn.getChrName().find("X") != string::npos) ||
			(toReturn.getChrName().find("Y") != string::npos) ){
			found=false;
		    }
		}
	    }
	}
    }

    return toReturn;
}
