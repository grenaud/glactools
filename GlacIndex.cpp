/*
 * GlacIndex
 * Date: Jul-25-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacIndex.h"



GlacIndex::GlacIndex(){

}

GlacIndex::~GlacIndex(){

}


string GlacIndex::usage() const{
        
    return string(string("")+"glactools index [options] <gl|ac file>\n"+
		  "\n"+
		  "\tThis program will index a <gl|ac> file given that is was bgzipped\n"
		  );

}

int GlacIndex::run(int argc, char *argv[]){
    string bgzf_file  = string(argv[argc-1]);
    BGZF * fpBGZF = bgzf_open(bgzf_file.c_str(), "r");

    if(fpBGZF == NULL){
	cout<<"ERROR "<< bgzf_file<<endl;
        exit (1);
    }

    int has_EOF = bgzf_check_EOF(fpBGZF);
    if (has_EOF == 0) {
        cerr<<"No EOF marker, likely due to an I/O error "<<endl;
    }else{
    	//cerr<<"EOF found"<<endl;
    }

    //detect BAM\1
    char bufferHeader[4];
    int magic_len = bgzf_read(fpBGZF, bufferHeader, 4);
    if( (magic_len != 4) 
        || 
	bufferHeader[0] != 'B' 
	||
	bufferHeader[1] != 'A' 
	||
	bufferHeader[2] != 'M' 
	||
	bufferHeader[3] != 1) {
        cerr<<"GlacIndex: Invalid magic number "<<magic_len<<"\t"<<bufferHeader<<endl;
        return 0;
    }else{
        //cerr<<"magic found"<<endl;
    }

    //ACF or GLF
    unsigned char  formattest [3];
    ssize_t     bytesread=0;
    bytesread = bgzf_read(fpBGZF, &formattest, 3);
    if(bytesread != 3){
	cerr<<"Error: GlacIndex tried to read "<< 3<<" bytes but got "<<bytesread<<endl;
	exit(1);
    }

    const string magicNumberACF="ACF";
    bool acfFound=true;
    const string magicNumberGLF="GLF";
    bool glfFound=true;
    
    for(unsigned int i=0;i<3;i++){
	if(formattest[i] != magicNumberACF[i]){
	    acfFound=false;
	}
    }
    
 	
    //Should find glf
    for(unsigned int i=0;i<3;i++){
	if(formattest[i] != magicNumberGLF[i]){
	    glfFound=false;
	}
    }
    
    if(!acfFound && !glfFound){
	cerr<<"Error: Did not find "<<magicNumberACF<<" nor "<<magicNumberGLF<<" found "<< formattest[0]<<formattest[1]<<formattest[2]<<endl;
	exit(1);
    }

    //size of fields in 1 byte
    char sizeBytes;
    bytesread = bgzf_read(fpBGZF, &sizeBytes, sizeof(sizeBytes));
    if(bytesread != sizeof(sizeBytes)){
	cerr<<"Error: GlacIndex tried to read "<< sizeof(sizeBytes) <<" bytes but got "<<bytesread<<endl;
	exit(1);
    }

    // \1
    char c;
    bytesread = bgzf_read(fpBGZF, &c, sizeof(c));
    if(bytesread != sizeof(c)){
	cerr<<"Error: GlacIndex tried to read "<< sizeof(c) <<" bytes but got "<<bytesread<<endl;
	exit(1);
    }
    
    
    uint32_t sizeHeader;
    size_t bytes = bgzf_read(fpBGZF, &sizeHeader, 4);
    if (bytes != 4){
	cerr<<"Could not read size of header"<<endl;
	return 0;
    }

    //cout<<"sizeHeader "<<sizeHeader<<endl;
    string header="";
    for(uint32_t i=0;i<sizeHeader;i++){
        char readChar;
	bytes = bgzf_read(fpBGZF, &readChar , 1);

	header += readChar;
	if (bytes != 1){
	    cerr<<"Could not read size header char #"<<i<<endl;
	    return 0;
	}
    }

    vector<string> allLinesHeader = allTokens(header,'\n');
    int32_t n_targets=0;
    for(unsigned l=0;l<allLinesHeader.size();l++){

	if(strBeginsWith(allLinesHeader[l],"#SQ\t")){
	    //cout<<allLinesHeader[l]<<endl;
	    n_targets++;
	}
    }

    uint32_t sizePops;
    bytes = bgzf_read(fpBGZF, &sizePops, 4);
    if (bytes != 4){
	cerr<<"Could not read size of population"<<endl;
	return 0;
    }else{
	//cout<<"found pop#"<<sizePops<<endl;
    }



    hts_idx_t *idx;

    //idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp), min_shift, n_lvls);
    //cerr<<"init hts_idx"<<endl;
    idx = hts_idx_init(n_targets, 1 , bgzf_tell(fpBGZF), 14, 5);
    
    //cout<<int(sizeBytes)<<endl;
    int sizeOneInd=0;
    if(acfFound){//ACf
	sizeOneInd = 2*int(sizeBytes)+1;
    }else{//GLf
	sizeOneInd = 3*int(sizeBytes)+1;
    }
    ssize_t sizeRecord = 8 + ( ssize_t(sizeOneInd) ) *(sizePops+2);

    char bufferRecord [sizeRecord];
    //cerr<<sizeRecord<<endl;
    //

    int ret;
    while( (bytesread = bgzf_read(fpBGZF, &bufferRecord, sizeRecord) ) ){
	if(bytesread==0)
	    break;

	if (bytesread != sizeRecord){
	    cerr<<"Could not read "<<sizeRecord<<" bytes, read "<<bytesread<<" instead"<<endl;
	    return 0;
	}

	uint16_t chr; //2bytes
        memcpy((char*)&chr,           bufferRecord+0,    sizeof(chr));//need correspondance to TID
        uint32_t coordinate ; //4bytes
        memcpy((char*)&coordinate,    bufferRecord+2,    sizeof(coordinate));

	ret = hts_idx_push(idx, chr, coordinate-1, coordinate, bgzf_tell(fpBGZF), 1);//last argument is format (fmt) 
	//cout<<chr<<"\t"<<coordinate<<endl;
	//exit(1);
	if(ret < 0){
	    cerr<<"GlacIndex: file "<<bgzf_file <<" is likely unsorted"<<endl;
	}
	// ret = hts_idx_push(idx, b->core.tid, b->core.pos, bam_endpos(b), bgzf_tell(fp), !(b->core.flag&BAM_FUNMAP));
        // if (ret < 0) goto err; // unsorted	
    }
    // char data;
    // int64_t tellPos_;

    hts_idx_finish(idx, bgzf_tell(fpBGZF));
    //cout<<"done"<<endl;
    bgzf_close(fpBGZF);//close file read

    
    //ret = hts_idx_save_as(idx, fn, fnidx, 1);//last argument is format (fmt) 
	    

    //write index
    //fp = bgzf_open(fnidx, (fmt == HTS_FMT_BAI)? "wu" : "w");
    string indexFile = bgzf_file+".bai";
    if(isFile(indexFile)){
	cerr<<"GlacIndex: WARNING overwriting "<<indexFile<<endl;
    }
    //BGZF fpIndex = bgzf_open(  , (fmt == HTS_FMT_BAI)? "wu" : "w"); 
    //cerr<<"writing index to indexFile..."<<indexFile<<endl;
    //BGZF * fpIndex = bgzf_open( indexFile.c_str() , "w");

    //bgzf_write(fpIndex, "BAI\1", 4);
    //check( hts_idx_save_core(idx, fpIndex, 1) );
    int returnCodeSave= hts_idx_save_as(idx, bgzf_file.c_str(), indexFile.c_str(), 1);
    //int returnCodeSave= hts_idx_save(idx, bgzf_file.c_str(),  1);
    if(returnCodeSave!=0){
	cerr<<"Problem writing index\n"<<endl;
	exit(1);
    }else{
	cerr<<"GlacIndex: index written to "<<indexFile<<endl;
    }
    //bgzf_close(fpIndex);

    hts_idx_destroy(idx);//destroy index
    
    
    return 0;


    return 0;
}


