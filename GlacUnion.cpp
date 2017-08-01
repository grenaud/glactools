/*
 * GlacUnion
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacUnion.h"

GlacUnion::GlacUnion(){

}

GlacUnion::~GlacUnion(){

}


string GlacUnion::usage() const{

    
    return string("glactools") +" union [OPTIONS] <acf|glf file1> <acf|glf file2> .. "+"\n"+
                  "\nThis program returns the intersection of ACF/GLF files given that they were from the same genome assembly\n"+    
		  "and prints to STDOUT. It will skip triallelic sites\n"+
	          "Files have to be from the same organism\n"+
                  "Options:\n"+
		  "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"
             	;


}



int GlacUnion::run(int argc, char *argv[]){

    
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
	cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;


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
    header<<"#GLACUNION:"<<endl;
    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	header<<"#UNIONFILE#"<<(i+1)<<endl;
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
    bool atLeastOneHasData=true;///should be true

    // vector<bool>  hasCoordinate (vectorOfGP.size(),true);//dummy value

   GlacWriter * gw=NULL;

   gw = new GlacWriter(numberOfPops,
		       isGL, //gp.isGLFormat(),
		       isGL?1:2,//gp.isACFormat()?2:1,
		       1,//compression threads
		       uncompressed);
   if(!gw->writeHeader(header.str())){
       cerr<<"GlacIntersect: error writing header "<<endl;
       exit(1);
   }
    

    bool stayLoop=true;


    while(stayLoop){
	if(!atLeastOneHasData ){
	    stayLoop=false;
	    break;
	}

#ifdef DEBUG
	cerr<<"coordCurrent "<<chr1<<"\t"<<coordCurrent<<endl;
#endif

	vector<bool> hasCoordinate (vectorOfGP.size(),false);
	bool atLeastOneHasCoordinate=false;
	for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	    if(hasData[i]){
		if(coordCurrent  == vecAlleleRecords[i]->coordinate){
		    hasCoordinate[i]=true;
		    atLeastOneHasCoordinate=true;
		}
	    }
	}
	// cout<<vectorToString(hasData,"-")<<endl;

	
	//we print
	if(atLeastOneHasCoordinate){
	    
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
	     atLeastOneHasData=false;

	     for(unsigned int i=0;i<vectorOfGP.size();i++){ 
		 
		 if(hasData[i] ){


		     //only get data from those with the coordinate
		     if(hasCoordinate[i]){
			 hasData[i]  =  vectorOfGP[i]->hasData();
			 if(hasData[i]){
			     atLeastOneHasData=true;
			     vecAlleleRecords[i] = vectorOfGP[i]->getData() ;
			 }
		     }else{
			 atLeastOneHasData=true;//still one with data
		     }

		}
	     }


	     // cout<<"seek "<<vectorToString(hasData,"-")<<"\t"<<atLeastOneHasData<<endl;

	     if(!atLeastOneHasData)
		 continue;

	     //coordCurrent=0;
	     //Try to find the record that is the most behind
	     bool needToSetCoord=true;
	     for(unsigned int i=0;i<vectorOfGP.size();i++){ 

		if(hasData[i]  ){		    
		    if(needToSetCoord){

			chr1          = vecAlleleRecords[i]->chri;
			coordCurrent  = vecAlleleRecords[i]->coordinate;

			needToSetCoord=false;
		    }else{
			if(chr1 == vecAlleleRecords[i]->chri){
			    coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
			}else{

			    int chrcmp = compare2ChrsU(chr1, vecAlleleRecords[i]->chri);
			    if(chrcmp == -1){// the current record is ahead, do nothing
				
				
			    }else{
				if(chrcmp== 1){ //the current record is behind chromosome-wise, reposition there.
				    chr1          = vecAlleleRecords[i]->chri;
				    coordCurrent  = vecAlleleRecords[i]->coordinate;

				}else{
				    cerr<<"Invalid state"<<endl;
				    return 1;  				    
				}
			    }
			}
		    }
		    // if(i==0){
		    // 	coordCurrent  = vecAlleleRecords[i]->coordinate;
		    // }else{
		    // 	coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
		    // }	
		}
	     }

	    continue;

	}else{//not at least one has coordinate
	     cerr<<"Invalid state"<<endl;
	     return 1;  
	}

	
    }


    for(unsigned int i=0;i<vectorOfGP.size();i++){ 
	delete(vectorOfGP[i]);
    }    


    cerr<<"Program glactools union terminated gracefully"<<endl;
    delete(gw);






    return 0;
}

