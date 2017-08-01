/*
 * GlacWriter
 * Date: Jul-27-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacWriter.h"

GlacWriter::GlacWriter(uint32_t sizePops_,bool glFormat_,int bytesForRecord_,int compressionThreads,bool uncompressed_){
    sizePops       = sizePops_;
    glFormat       = glFormat_;
    bytesForRecord = bytesForRecord_;
    uncompressed   = uncompressed_;

    if(glFormat){
	if(bytesForRecord!=1){
	    cerr<<"For GLF, the number of bytes should be 1"<<endl;
	    exit(1);
	}
    }else{
	if(bytesForRecord!=2 && 
	   bytesForRecord!=3){
	    cerr<<"For ACF, the number of bytes should be 2 or 3"<<endl;
	    exit(1);
	}
    }
    //cout<<bytesForRecord<<endl;

    string bgzf_file = "/dev/stdout";
    fpBGZF    = NULL;

    if(uncompressed){

    }else{//compressed
	fpBGZF = bgzf_open(bgzf_file.c_str(), "w");
	//fpBGZF = bgzf_dopen(bgzf_file.c_str(), "w");

	if(compressionThreads>1)
	    bgzf_mt(fpBGZF, compressionThreads, 256);

	if (fpBGZF == NULL) { // region invalid or reference name not found
	    cerr<<"Cannot write to "<<bgzf_file<<endl;
	    exit(1);
	}else{		
	}
    }

    
}

GlacWriter::~GlacWriter(){
    if(!uncompressed){
	if(bgzf_close(fpBGZF) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }  
    }

}

bool GlacWriter::writeHeader(string const & header) const{

    char bammagicstr [4] = {'B','A','M','\1'};    
    //cerr<<uncompressed<<endl;
    if(uncompressed){
	if(write(1,&bammagicstr,sizeof(bammagicstr)) == -1 ){   cerr<<"Write error"<<endl;            return false;   }     
    }else{
	if( bgzf_write(fpBGZF, &bammagicstr,sizeof(bammagicstr)) != sizeof(bammagicstr) ){   cerr<<"Write error"<<endl;   return false;   }     
    }
    
    if(!glFormat){
	char magicstr [5] = {'A','C','F', char(bytesForRecord) ,'\1'};
	if(uncompressed){
	    if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return false;   }     
	}else{
		    if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;  return false;   }     
	}
    }else{		
	char magicstr [5] = {'G','L','F', char(bytesForRecord) ,'\1'};   
	if(uncompressed){
	    if(write(1,&magicstr,sizeof(magicstr)) == -1 ){   cerr<<"Write error"<<endl;            return false;   }     
	}else{
	    if( bgzf_write(fpBGZF, &magicstr,sizeof(magicstr)) != sizeof(magicstr) ){   cerr<<"Write error"<<endl;  return false;   } 
		}
    }
    
    uint32_t sizeHeader = header.size();
    if(uncompressed){
	if(write(1,&sizeHeader,sizeof(sizeHeader)) == -1 ){   cerr<<"Write error"<<endl;        return false;   } 
    }else{
	if( bgzf_write(fpBGZF, &sizeHeader,sizeof(sizeHeader)) != sizeof(sizeHeader) ){   cerr<<"Write error"<<endl; return false;   }     
    }
    
    for(uint32_t i=0;i<sizeHeader;i++){
	char towrite=char(header[i]);
	if(uncompressed){
	    if(write(1,&towrite,sizeof(towrite)) == -1 ){   cerr<<"Write error"<<endl;  return false;   } 
	}else{
	    if( bgzf_write(fpBGZF, &towrite,sizeof(towrite)) != sizeof(towrite) ){   cerr<<"Write error"<<endl; return false;   }     
	}
    }
    
    //uint32_t sizePops=gp.getSizePops();
    if(uncompressed){
	if(write(1,&sizePops,sizeof(sizePops)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
    }else{
	if( bgzf_write(fpBGZF, &sizePops,sizeof(sizePops)) != sizeof(sizePops) ){   cerr<<"Write error"<<endl;            return false;   }     
    }

    return true;
}

bool GlacWriter::writeAlleleRecord(const AlleleRecords * toWrite) const{
    

    if(uncompressed){
	if(      write(1,      &toWrite->chri,sizeof(toWrite->chri)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
    }else{
	if( bgzf_write(fpBGZF, &toWrite->chri,sizeof(toWrite->chri)) != sizeof(toWrite->chri) ){   cerr<<"Write error"<<endl;   return false;   }  
    }

#ifdef DEBUGWRITEAR
    cerr<<"chri "<<toWrite->chri<<endl;    
#endif

    //write coordinate 4 bytes
    //uint32_t c = uint32_t(toprint->getPosition());
    //if(write(1,&c,    sizeof(c))  == -1 )    {   cerr<<"Write error"<<endl;         return false;   }
    if(uncompressed){
	if(       write(1,     &toWrite->coordinate,sizeof(toWrite->coordinate)) == -1 ){   cerr<<"Write error"<<endl;      return false;   } 
    }else{
	if( bgzf_write(fpBGZF, &toWrite->coordinate,sizeof(toWrite->coordinate)) != sizeof(toWrite->coordinate) ){   cerr<<"Write error"<<endl;  return false;   }  
    }
#ifdef DEBUGWRITEAR
    cerr<<"coordinate "<<toWrite->coordinate<<endl;    		
#endif

    uint8_t  tempCh;
    uint16_t tempSh;

#ifdef DEBUGWRITEREFALT
    //cerr<<"ref "<<"NACGT"[tempCh]<<endl;
    cerr<<"ref "<<toWrite->ref<<"\t"<<int(tempCh)<<endl;
#endif
    
    //ref 1 byte
    tempCh= uint8_t(base2int(toWrite->ref));
    tempCh=tempCh<<4;//shift by 4 bits
    //if(write(1,&tempCh,sizeof(tempCh))  == -1 ){   cerr<<"Write error"<<endl;           return false;   }
    


    // if(uncompressed){
    // 	if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
    // }else{
    // 	if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;    return false;   }  
    // }
    
    //alt 1 byte
    //tempCh= uint8_t(base2int(toWrite->alt));
#ifdef DEBUGWRITEREFALT
    //cerr<<"alt "<<"NACGT"[tempCh]<<endl;
    cerr<<"ref+alt "<<toWrite->alt<<"\t"<<int(tempCh)<<endl;
#endif

    tempCh = tempCh | uint8_t(base2int(toWrite->alt));


    if(uncompressed){
	if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
    }else{
	if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;    return false;   }  
    }
    

    // cout<<"sizePops "<<sizePops<<endl;
    // exit(1);
#ifdef DEBUGWRITEAR
    cerr<<"glformat "<<glFormat<<endl;
#endif
    
    if(!glFormat){//not GLF so ACF
	if(int(toWrite->vectorAlleles->size()) != int(sizePops+2)){
	    cerr<<"Discrepancy between vector of alleles "<<int(toWrite->vectorAlleles->size())<<" and given number of pops "<< int(sizePops+2)<<endl;
	    return false;
	}

#ifdef DEBUGWRITEAR
	cerr<<"not gl"<<endl;
#endif
	if(bytesForRecord==2){

	    uint8_t  toWriteCpG=0; 
	    
	    for(uint32_t j=0;j<(sizePops+2);j++){//plus 2 for root anc
		//2 bytes
		tempSh= uint16_t( toWrite->vectorAlleles->at(j).getRefCount() );
#ifdef DEBUGWRITEAR
		cerr<<j<<" ref count"<<tempSh<<endl;
#endif

		if(uncompressed){
		    if(write(1,&tempSh,sizeof(tempSh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempSh,sizeof(tempSh)) != sizeof(tempSh) ){   cerr<<"Write error"<<endl;            return false;   }  
		}

		//2 bytes
		tempSh= uint16_t( toWrite->vectorAlleles->at(j).getAltCount() );
#ifdef DEBUGWRITEAR
		cerr<<j<<" alt count"<<tempSh<<endl;
#endif

		
		if(uncompressed){
		    if(write(1,&tempSh,sizeof(tempSh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempSh,sizeof(tempSh)) != sizeof(tempSh) ){   cerr<<"Write error"<<endl;            return false;   }  
		}

		//1 byte
		tempCh= uint8_t( toWrite->vectorAlleles->at(j).getIsCpg() );

#ifdef DEBUGWRITEAR
		cerr<<j<<" alt count"<<tempCh<<endl;
#endif
		
		if(uncompressed){
		    if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
		}else{
		    if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return false;   }  
		}

	    
	    }//end for each pop
	}else{//bytes for record=2
	    cerr<<"not yet implemented"<<endl;
	    exit(1);
	}
    }//end acformat
    else{

	if(int(toWrite->vectorGLs->size()) != int(sizePops+2)){
	    cerr<<"Discrepancy between vector of alleles "<<int(toWrite->vectorAlleles->size())<<" and given number of pops "<< int(sizePops+2)<<endl;
	    return false;
	}

	for(uint32_t j=0;j<(sizePops+2);j++){//root anc
	    //1 byte
	    tempCh= toWrite->vectorGLs->at(j).getrrGL() ;
	    if(uncompressed){
		if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
	    }else{
		if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return false;   }  
	    }

	    //1 byte
	    tempCh= toWrite->vectorGLs->at(j).getraGL() ;
	    if(uncompressed){
		if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
	    }else{
		if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return false;   }  
	    }
	    
	    //1 byte
	    tempCh= toWrite->vectorGLs->at(j).getaaGL();
	    if(uncompressed){
		if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
	    }else{
		if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return false;   }  
	    }
	    
	    tempCh= uint8_t( toWrite->vectorGLs->at(j).getIsCpg() );
	    if(uncompressed){
		if(write(1,&tempCh,sizeof(tempCh)) == -1 ){   cerr<<"Write error"<<endl;          return false;   } 
	    }else{
		if( bgzf_write(fpBGZF, &tempCh,sizeof(tempCh)) != sizeof(tempCh) ){   cerr<<"Write error"<<endl;            return false;   }  
	    }
	}//end for each pop
	
    }//end glformat
    return true;
}//end write
