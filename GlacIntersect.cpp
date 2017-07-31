/*
 * GlacIntersect
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacIntersect.h"

GlacIntersect::GlacIntersect(){

}

GlacIntersect::~GlacIntersect(){

}


string GlacIntersect::usage() const{

    
    return string("glactools") +" glacintersect [OPTIONS] <acf|glf file1> <acf|glf file2> .. "+"\n"+
                  "\nThis program returns the intersection of ACF/GLF files given that they were from the same genome assembly\n"+    
		  "and prints to STDOUT. It will skip triallelic sites\n"+
	          "Files have to be from the same organism\n"+
                  "Options:\n"+
		  "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"
             	;

}



int GlacIntersect::run(int argc, char *argv[]){

    bool force=false;    


    if(argc < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr<<"usage: "<<usage()<<endl;
	return 1;       
    }
   int lastOpt=1;
			      

    //last arg is program name
    for(int i=1;i<(argc-1);i++){ 
	if((string(argv[i]) == "-")  ){
            lastOpt=i;
            break;	    
	}

        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }                               

	if(strcmp(argv[i],"-f") == 0 ){
	    force=true;
            continue;
	}

        if(string(argv[i]) == "-u"){
	    uncompressed=true;
            continue;
        }
    }

    string defLine;
    string sqLines;
    bool isGL=false;
    vector<GlacParser * > vectorOfGP;
    int numberOfPops=0;
    for(int i=lastOpt;i<(argc);i++){ 

	GlacParser * gp = new GlacParser(string(argv[i]));

	numberOfPops += int(gp->getSizePops());

	if(i==lastOpt){
	    sqLines = gp->getHeaderSQ();
	    isGL = gp->isGLFormat();

	}else{
	    if(isGL != gp->isGLFormat()){
		cerr<<"GlacIntersect: Error file "<<string(argv[i])<<" is not not in the same format as "<<string(argv[lastOpt])<<endl;
		return 1;
	    }

	    if(    sqLines != gp->getHeaderSQ() ){
		cerr<<"GlacIntersect: Error file "<<string(argv[i])<<" does not have the same SQ lines as "<<string(argv[lastOpt])<<endl;
		return 1;
	    }
	}
	
	vectorOfGP.push_back(gp);
    }    



    stringstream header;
    if(isGL)
	header<<"#GLF"<<endl;    	
    else
	header<<"#ACF"<<endl;    	

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    header<<"#PG:"<<programLine<<endl;
    header<<"#GITVERSION: "<<returnGitHubVersion(argv[-1],"")<<endl;
    header<<"#DATE: "<<getDateString()<<endl;
    header<<"#GLACINTERSECT:"<<endl;
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	header<<"#INTERSECTFILE#"<<(i+1)<<endl;
	header<<""<<vectorOfGP[i]->getHeaderNoSQNoDefline("#\t")<<endl;
    }
    header<<sqLines<<endl;
    //header<<
    //    bool atLeastOneHasData;
    vector<bool> hasData;
    vector<int> popSizePerFile;
    vector<AlleleRecords *> vecAlleleRecords;
    //string chr1;
    uint16_t chr1;
    unsigned int coordCurrent;


    
    string defline=initFiles(vectorOfGP,
			     // atLeastOneHasData,
			     hasData,
			     popSizePerFile,
			     vecAlleleRecords,
			     chr1,
			     coordCurrent);
    header<<defline<<endl;
    vector<bool>  hasCoordinate (vectorOfGP.size(),true);//dummy value

   GlacWriter * gw=NULL;

   gw = new GlacWriter(numberOfPops,
		       isGL, //gp.isGLFormat(),
		       isGL?1:2,//gp.isACFormat()?2:1,
                        uncompressed);
   if(!gw->writeHeader(header.str())){
       cerr<<"GlacIntersect: error writing header "<<endl;
       exit(1);
   }
    

    bool stayLoop=true;


    while(stayLoop){

#ifdef DEBUG
	cerr<<"coordCurrent "<<chr1<<":"<<coordCurrent<<endl;
#endif

	
	bool allHaveCoordinate=true;

	for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	    if(hasData[i]){		
		if( (chr1          != vecAlleleRecords[i]->chri      ) ||
		    (coordCurrent  != vecAlleleRecords[i]->coordinate) ){
		    allHaveCoordinate=false;
		}
	    }else{
		stayLoop=false;
		break;
	    }
	}

	
	//we print
	if(allHaveCoordinate){
#ifdef DEBUG
	cerr<<"same coordCurrent "<<chr1<<":"<<coordCurrent<<endl;
#endif

	    
	    printAllele(vectorOfGP,
			hasData,
			hasCoordinate,
			popSizePerFile,
			vecAlleleRecords,
			chr1,
			coordCurrent,
			gw,
			isGL,
			force);


	    // 	seekdata:
	    allHaveCoordinate=false;
	    
	    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
		 
		if(!hasData[i] ){
		    cerr<<"Invalid state"<<endl;
		    return 1;
		}
		hasData[i]  =  vectorOfGP[i]->hasData();
		if(hasData[i]){
		    vecAlleleRecords[i] = vectorOfGP[i]->getData() ;
		}else{
		    stayLoop=false;
		    break;
		}
		
	    }

	    
	    //all have had getData called, we need to reposition to the maximum coord
	    bool needToSetCoord=true;
	    for(unsigned int i=0;i<vectorOfGP.size();i++){ 

#ifdef DEBUG
		cerr<<needToSetCoord<<"\t"<<i<<"\t"<<chr1<<":"<<coordCurrent<<endl;
#endif

		 if(needToSetCoord){
		     chr1          = vecAlleleRecords[i]->chri;
		     coordCurrent  = vecAlleleRecords[i]->coordinate;
		     needToSetCoord=false;
		 }else{
		     int chrcmp = compare2ChrsU(chr1, vecAlleleRecords[i]->chri);

		     if(chrcmp != 0){//chromosomes are not equal

			 if(chrcmp == -1){ //chr1 < vecAlleleRecords[i]->chr
			     chr1          = vecAlleleRecords[i]->chri;
			     coordCurrent  = vecAlleleRecords[i]->coordinate;			     
			 }

			 if(chrcmp ==  1){ // chr1 > vecAlleleRecords[i]->chr
			      // chr1          = chr1;
			      // coordCurrent  = vecAlleleRecords[i]->coordinate;
			 }

			 

		     }else{
			 coordCurrent  = max(coordCurrent,vecAlleleRecords[i]->coordinate); //chromosomes are equal, jump to the max
		     }


		     // if( chr1          != vecAlleleRecords[i]->chr ){//need to skip to a new chr
		     // 	 if(coordCurrent  > vecAlleleRecords[i]->coordinate){//the diff chr is probably a new chr
		     // 	     chr1          = chr1;
		     // 	     coordCurrent  = vecAlleleRecords[i]->coordinate;
		     // 	 }
		     // }else{
		     // 	 if( chr1          == vecAlleleRecords[i]->chr ){
		     // 	     coordCurrent  = max(coordCurrent,vecAlleleRecords[i]->coordinate);
		     // 	 }else{//chr1 is greater than vecAlleleRecords[i]
		     // 	     //	we will reposition vecAlleleRecords[i] there 
		     // 	 }
			 
		     // }
		 }

	     }

	    continue;

	}else{
	     // cerr<<"Invalid state"<<endl;
	     // return 1;  
	    

#ifdef DEBUG
	    cerr<<"current "<<chr1<<"\t"<<coordCurrent<<endl;
#endif
	    
	    for(unsigned int i=0;i<vectorOfGP.size();i++){ 

#ifdef DEBUG
		cerr<<"coord["<<i<<"] "<<vecAlleleRecords[i]->chri<<"\t"<<vecAlleleRecords[i]->coordinate<<endl;
#endif

		//record [i] is ahead, re-position there
		if( ((chr1          == vecAlleleRecords[i]->chri) && 
		     (coordCurrent  < vecAlleleRecords[i]->coordinate)) ){ //overshot [i], repositioning there
		    chr1          = vecAlleleRecords[i]->chri       ;
		    coordCurrent  = vecAlleleRecords[i]->coordinate;
		    continue;
		}

		int chrcmp = compare2ChrsU(chr1, vecAlleleRecords[i]->chri);

		//record [i] is ahead, re-position there
		if( chrcmp  ==  -1    ){// chr1 < vecAlleleRecords[i]->chr
		//if( chr1          != vecAlleleRecords[i]->chr         ){// chr1 < vecAlleleRecords[i]->chr

		    //if(coordCurrent  > vecAlleleRecords[i]->coordinate)  { //the different chr is probably a new chr
		    chr1          = vecAlleleRecords[i]->chri      ;
		    coordCurrent  = vecAlleleRecords[i]->coordinate;
		    continue;
		    //}

		}

		//if( (chr1         == vecAlleleRecords[i]->chr) &&
		//in sync
		if( (chrcmp        == 0) &&
		    (coordCurrent == vecAlleleRecords[i]->coordinate)){ //fine
		}
		    
		//if( ((chr1          != vecAlleleRecords[i]->chr) && (coordCurrent <  vecAlleleRecords[i]->coordinate)) 
		//record [i] is running behind 
		if( (chrcmp == 1)  // chr1 > vecAlleleRecords[i]->chr
		    ||
		    ((chr1          == vecAlleleRecords[i]->chri) && (coordCurrent >  vecAlleleRecords[i]->coordinate)) ){ 
		    hasData[i]  =  vectorOfGP[i]->hasData();
		    if(hasData[i]){
			vecAlleleRecords[i] = vectorOfGP[i]->getData() ;
		    }else{
			stayLoop=false;
			break;
		    }
		}
		
	    }

	    
	}//end different coord

	
    }//end main loop

    // finish:
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 

	delete(vectorOfGP[i]);
    }    
    delete(gw);

    cerr<<"Program glactools intersect terminated gracefully"<<endl;

    return 0;


    return 0;
}

