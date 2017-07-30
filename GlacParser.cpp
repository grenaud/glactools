/*
 * GlacParser
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacParser.h"



static int glac_readrecACF2b(BGZF *fp, void *ignored, void *bv, int *tid, int *beg, int *end){
    // bam1_t *b = bv;
    // int ret;
    AlleleRecords * ptrar = (AlleleRecords *) (bv);

    uint32_t sizePops=ptrar->getSizePops();
    //                    5b record (2b+2b+1)
    size_t sizeSingleAC = (2*sizeBytesACF+1);
    //                  8b base,5b record (2b+2b+1)
    size_t sizeRecord = 8+ (sizeSingleAC)*(sizePops+2);//+2 for root/anc
    char buffer [sizeRecord];
    // cout<<"sp "<<sizePops<<"\t"<<sizeRecord<<endl;
    // exit(1);


    //cout<<"glac_readrec  ignored "<<ignored<<" bv "<<bv<<" tid "<<*tid<<" beg "<<*beg<<" end "<<*end<<endl;
    //printf("glac_readrec %d\n",fp);

    //cout<<ptrar<<"\t"<<sizeof(ptrar->chri)<<endl;
    //myFilezipped.read((char*)&buffer,        sizeof(buffer));
    //COPY data
    //cout<<sizeRecord<<endl;
    size_t bytesread = bgzf_read(fp, &buffer, sizeRecord);
    if(bytesread != sizeRecord){
	cerr<<"GlacParser: tried to read "<<sizeRecord<<" but got "<<bytesread<<endl;
	exit(1);
    }
    //cout<<"bytesread "<<bytesread<<endl;        
    // if ((ret = bam_read1(fp, b)) >= 0) {
    //     *tid = b->core.tid;
    //     *beg = b->core.pos;
    //     *end = bam_endpos(b);
    // }
    size_t offset = 0;

    //uint16_t chr; //2
    //memcpy((char*)&chr,           buffer+offset,    sizeof(chr));
    memcpy((char*)&ptrar->chri,        buffer+offset,    sizeof(ptrar->chri));
    //cout<<"chr "<<ptrar->chri<<endl;


    offset+=sizeof(ptrar->chri);
    //uint32_t coordinate ; //4
    //memcpy((char*)&coordinate,    buffer+offset,    sizeof(coordinate));
    memcpy((char*)&ptrar->coordinate,  buffer+offset,    sizeof(ptrar->coordinate));
    //cout<<"coor "<<ptrar->coordinate<<endl;

    offset+=sizeof(ptrar->coordinate);

    char ref; //1
    memcpy((char*)&ref,           buffer+offset,    sizeof(ref));
    //memcpy((char*)&ptrar->ref,         buffer+offset,    sizeof(ptrar->ref));
    ptrar->ref = "NACGT"[(unsigned char)ref];
    offset+=sizeof(ptrar->ref);
    //cout<<"ref  "<<"NACGT"[ptrar->ref]<<endl;
    char alt; //1
    memcpy((char*)&alt,           buffer+offset,    sizeof(alt));
    //memcpy((char*)&ptrar->alt,         buffer+offset,    sizeof(ptrar->alt));
    ptrar->alt = "NACGT"[(unsigned char)alt];
    offset+=sizeof(ptrar->alt);    
    // cout<<"alt  "<<"NACGT"[ptrar->alt]<<endl;
    // cout<<"chr "<<ptrar->chri<<endl;    
    //    cout<<"chr "<<ptrar->chri<<"\t";
    //cout<<ptrar->coordinate<<endl;
    //cout<<"NACGT"[ptrar->ref]<<","<<"NACGT"[ptrar->alt]<<"\t";

    *beg=int(ptrar->coordinate);
    *end=int(ptrar->coordinate);
    *tid=int(ptrar->chri);

    ptrar->vectorAlleles = new vector<SingleAllele>();
    
    for(unsigned j=0;j<(sizePops+2);j++){
	//cout<<j<<" "<<sizePops<<endl;
	uint16_t refC; //2
	uint16_t altC; //2
	char     cpgC; //1
	//5 = 2+2+1
	memcpy((char*)&refC,       buffer+8 +5*j, sizeof(refC));
	memcpy((char*)&altC,       buffer+10+5*j, sizeof(altC));
	memcpy((char*)&cpgC,       buffer+12+5*j, sizeof(cpgC));

	// cout<<refC<<endl;
	// cout<<altC<<endl;
	// cout<<(cpgC==1)<<endl;
	 // if(j < (sizePops -1 ))
	 //     cout<<"\t";

	// memcpy((char*)&ptrar->vectorAlleles->[j=refC,       buffer+8 +5*j, sizeof(ptrar->refC));
	// memcpy((char*)&ptrar->altC,       buffer+10+5*j, sizeof(ptrar->altC));
	// memcpy((char*)&ptrar->cpgC,       buffer+12+5*j, sizeof(ptrar->cpgC));
	//ptrar->vectorAlleles->at(j).setRefCount(refC);
	//ptrar->vectorAlleles->at(j).setAltCount(altC);
	//ptrar->vectorAlleles->at(j).setIsCpg(   cpgC);
	SingleAllele sa (int( refC ),
			 int( altC ),
			 (cpgC == 1) );
	//cout<<j<<"\t"<<sa<<endl;
	ptrar->vectorAlleles->push_back(sa);		         
    }
    //cout<<endl;
    //cout<<"test#"<<*ptrar<<"#"<<endl;
    

    return sizeRecord;
}






static int glac_readrecGLF1b(BGZF *fp, void *ignored, void *bv, int *tid, int *beg, int *end){
    // bam1_t *b = bv;
    // int ret;
    AlleleRecords * ptrar = (AlleleRecords *) (bv);

    uint32_t sizePops=ptrar->getSizePops();
    //                     4b record (1+1+1+1)
    size_t sizeSingleGL = (3*sizeBytesGLF+1);
    //                  8b base,4b record (1+1+1+1)
    size_t sizeRecord = 8+ sizeSingleGL*(sizePops+2);//+2 for root/anc
    char buffer [sizeRecord];


    //cout<<"glac_readrec  ignored "<<ignored<<" bv "<<bv<<" tid "<<*tid<<" beg "<<*beg<<" end "<<*end<<endl;
    //printf("glac_readrec %d\n",fp);

    //cout<<ptrar<<"\t"<<sizeof(ptrar->chri)<<endl;
    //myFilezipped.read((char*)&buffer,        sizeof(buffer));
    //COPY data
    //cout<<sizeRecord<<endl;
    size_t bytesread = bgzf_read(fp, &buffer, sizeRecord);
    if(bytesread != sizeRecord){
	cerr<<"GlacParser: tried to read "<<sizeRecord<<" but got "<<bytesread<<endl;
	exit(1);
    }
    //cout<<"bytesread "<<bytesread<<endl;        
    // if ((ret = bam_read1(fp, b)) >= 0) {
    //     *tid = b->core.tid;
    //     *beg = b->core.pos;
    //     *end = bam_endpos(b);
    // }
    size_t offset = 0;

    //uint16_t chr; //2
    //memcpy((char*)&chr,           buffer+offset,    sizeof(chr));
    memcpy((char*)&ptrar->chri,        buffer+offset,    sizeof(ptrar->chri));
    //cout<<"chr "<<ptrar->chri<<endl;


    offset+=sizeof(ptrar->chri);
    //uint32_t coordinate ; //4
    //memcpy((char*)&coordinate,    buffer+offset,    sizeof(coordinate));
    memcpy((char*)&ptrar->coordinate,  buffer+offset,    sizeof(ptrar->coordinate));
    //cout<<"coor "<<ptrar->coordinate<<endl;

    offset+=sizeof(ptrar->coordinate);

    char ref; //1
    memcpy((char*)&ref,           buffer+offset,    sizeof(ref));
    //memcpy((char*)&ptrar->ref,         buffer+offset,    sizeof(ptrar->ref));
    ptrar->ref = "NACGT"[(unsigned char)ref];
    offset+=sizeof(ptrar->ref);
    //cout<<"ref  "<<"NACGT"[ptrar->ref]<<endl;
    char alt; //1
    memcpy((char*)&alt,           buffer+offset,    sizeof(alt));
    //memcpy((char*)&ptrar->alt,         buffer+offset,    sizeof(ptrar->alt));
    ptrar->alt = "NACGT"[(unsigned char)alt];
    offset+=sizeof(ptrar->alt);    
    // cout<<"alt  "<<"NACGT"[ptrar->alt]<<endl;
    // cout<<"chr "<<ptrar->chri<<endl;    
    //    cout<<"chr "<<ptrar->chri<<"\t";
    //cout<<ptrar->coordinate<<endl;
    //cout<<"NACGT"[ptrar->ref]<<","<<"NACGT"[ptrar->alt]<<"\t";

    *beg=int(ptrar->coordinate);
    *end=int(ptrar->coordinate);
    *tid=int(ptrar->chri);

    ptrar->vectorGLs = new vector<SingleGL>();
    
    for(unsigned j=0;j<(sizePops+2);j++){
	//cout<<j<<" "<<sizePops<<endl;
	uint8_t rrC;//1b
	uint8_t raC;//1b
	uint8_t aaC;//1b
	char   cpgC; //1

	memcpy((char*)&rrC,        buffer+ 8 +4*j, sizeof(rrC));
	memcpy((char*)&raC,        buffer+ 9 +4*j, sizeof(raC));
	memcpy((char*)&aaC,        buffer+10 +4*j, sizeof(aaC));
	memcpy((char*)&cpgC,       buffer+11 +4*j, sizeof(cpgC));

	// cout<<refC<<endl;
	// cout<<altC<<endl;
	// cout<<(cpgC==1)<<endl;
	 // if(j < (sizePops -1 ))
	 //     cout<<"\t";

	// memcpy((char*)&ptrar->vectorAlleles->[j=refC,       buffer+8 +5*j, sizeof(ptrar->refC));
	// memcpy((char*)&ptrar->altC,       buffer+10+5*j, sizeof(ptrar->altC));
	// memcpy((char*)&ptrar->cpgC,       buffer+12+5*j, sizeof(ptrar->cpgC));
	//ptrar->vectorAlleles->at(j).setRefCount(refC);
	//ptrar->vectorAlleles->at(j).setAltCount(altC);
	//ptrar->vectorAlleles->at(j).setIsCpg(   cpgC);

	SingleGL gl (rrC,
		     raC,
		     aaC,			     
		     (cpgC == 1) );
	ptrar->vectorGLs->push_back(gl);

    }
    //cout<<endl;
    //cout<<"test#"<<*ptrar<<"#"<<endl;
    

    return sizeRecord;
}



GlacParser::GlacParser(string bgzf_file,string bgzf_fileidx,string chrName,int start_,int end_,bool justChr){
    //cout<<chrName<<"\t"<<start<<"\t"<<end<<endl;
    //for rand() ins SingleAllele
    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );


    

    //     rt =new ReadTabix (file,indexForFile,chrName,start,end);
    //     string headertemp=rt->getHeader();
    
    //     istringstream is(headertemp);
    bool isbgzip=bgzf_is_bgzf(bgzf_file.c_str());

    myFilezipped = bgzf_open(bgzf_file.c_str(), "r");

    //reading the header
    numberPopulations=0;
    populationNames=new vector<string>();

    parseHeader(myFilezipped);

    if(defline.empty()){
	cerr << "Error: GlacParser cannot get definition line" << bgzf_file <<endl;
	exit(1);
    }

    //loading index
    htsFile *fp = hts_open(bgzf_file.c_str(),"r");
    

    hts_idx_t *idx = sam_index_load(fp, bgzf_fileidx.c_str()); // load index
    if (idx == 0) { // index is unavailable
    	cerr<<"Cannot load index "<<bgzf_fileidx<<endl;
    	exit(1);
    }else{
    	//cerr<<"index loaded succesfully\n"; //need to have bai index
    }


    int tid=0;
    for(int i=0;i<int(chrKnown.size());i++){
	if(chrKnown[i] == chrName){
	    tid=i;
	    break;
	}
    }
     

    if(start_>end_){
	cerr<<"GlacParser: start coordinate "<<start_<<" is greater than "<<end_<<endl;
	exit(1);
    }

    int beg=(start_-1);//-1 to make it include the start coord
    int end=end_;
    if(justChr){
	beg=0;
	end=INT_MAX;
    }
    //hts_itr_t *iter;
    if(acFormat && sizeBytesACF==2)
	iter= hts_itr_query(idx, tid, beg, end, glac_readrecACF2b); 

    if(glFormat && sizeBytesGLF==1)
	iter= hts_itr_query(idx, tid, beg, end, glac_readrecGLF1b); //TODO change for glf

    //cout<<"TID "<<tid<<" beg "<<beg<<" end "<<end<<endl;

    if (iter == NULL) { // region invalid or reference name not found
    	cerr<<"GlacParser: query failed to "<<"TID "<<tid<<" beg "<<beg<<" end "<<end<<endl;
    	exit(1);
    }else{
    	//cerr<<"query success!"<<iter->finished<<endl;
    }

    if(bgzf_close(myFilezipped) != 0){
	cerr<<"GlacParser: Cannot close input file "<<endl;	    
    }

    //should we close htsfile?    
    myFilezipped = bgzf_open(bgzf_file.c_str(), "r");
    //int result;

    //exit(1);

    //cout<<defline<<endl;
    //     parseHeader(is);    
    //     if(defline.empty()){
    // 	cerr << "Error: GlacParser cannot get definition line"  <<endl;
    // 	exit(1);
    //     }
    
    numberOfTimesHasDataWasCalled=-1;
    binMode    = false;
    tabixMode  = true;
    textMode   = false;
    stringMode = false;
 }

//TODO to implement for parallele
// GlacParser::GlacParser(const vector<string> * dataToRead_,const vector<string> & populationNames_){
//     //populationNames = populationNames_;
//     populationNames = new vector<string>( populationNames_);
//     numberPopulations = populationNames->size();
//     dataToRead      = dataToRead_;
//     header="";
//     headerNoDefline="";

//     numbernew=0;
//     numberdel=0;
//     dataToReadInd=0;
//     numberOfTimesHasDataWasCalled=-1;

//     stringMode = true;
//     tabixMode  = false;
//     textMode   = false;
// }

GlacParser::GlacParser(string filename){
    header="";
    headerNoDefline="";

    defline="";


    numbernew=0;
    numberdel=0;

    bool isbgzip;

    //cout<<"parser:"<<filename<<"#"<<endl;
    //filename="-";
    bool openSTDIN = (filename == "-" || filename == "/dev/stdin");//is there a better way to do this?
    //cerr<<openSTDIN<<endl;
    //exit(1);
    if(openSTDIN){
	//cerr<<"stdin"<<endl;
	myFilezipped=bgzf_dopen(0, "r");
    }else{
	isbgzip     =bgzf_is_bgzf(filename.c_str());
	myFilezipped=bgzf_open(filename.c_str(), "r");
    }    

    if(myFilezipped == 0){
        cerr<<"Error: GlacParser failed to open file "<< filename <<endl;
        exit (1);
    }

    if(!openSTDIN && isbgzip){
	int has_EOF = bgzf_check_EOF(myFilezipped);
	if (has_EOF == 0) {
	    cerr<<"Warning: No EOF marker, likely due to an I/O error "<<endl;	
	}else{
	    //cerr<<"EOF found"<<endl;
	}    
    }
    // myFilezipped=new igzstream();
    // myFilezipped->open( filename.c_str() );


    //for rand() ins SingleAllele
    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );


    // if ( ! myFilezipped->good()) {
    // 	cerr << "Error: GlacParser failed to open file=" << filename <<endl;
    // 	exit(1);;
    // }

    //reading the header
    numberPopulations=0;
    populationNames=new vector<string>();

    parseHeader(myFilezipped);
    //exit(1);
    if(defline.empty()){
	cerr << "Error: GlacParser cannot get definition line" << filename <<endl;
	exit(1);
    }

    numberOfTimesHasDataWasCalled=-1;
    
    // tabixMode  = false;
    // textMode   = true;
    // stringMode = false;
}

GlacParser::~GlacParser(){
    //cerr<<"destructor GlacParser "<<tabixMode<<" "<<binMode<<" "<<textMode<<" "<<numberOfTimesHasDataWasCalled<<endl;

    if(numberOfTimesHasDataWasCalled == 0){//last called was getData
	//delete(allRecToReturn->vectorAlleles);
	delete(allRecToReturn);
	numberdel++;
    }

    //TODO
    if(tabixMode){
	// 	delete rt; //calling the destructor
     }

    if(binMode){
	if(bgzf_close(myFilezipped) != 0){
	    cerr<<"GlacParser: Cannot close input file "<<endl;	    
	}
    }

    delete(populationNames);
    if(textMode){
	delete(myFilezipped);
    }

    // else
    // 	delete(myFile);
}

void GlacParser::parseHeader(BGZF *bg){
    //bool firstLine=true;
    //string line;
    //cout<<"parsed "<<textMode<<endl;
    //if(textMode){
    //cout<<"parsed"<<endl;
    const string magicNumberBAM="BAM\1";
    const string magicNumberACFt="#ACF";
    const string magicNumberGLFt="#GLF";
    string firstLineHeadertxtmode="";
    unsigned char  bamtest [4];
    const ssize_t     bamtestlength  = 4;
    ssize_t bytesread;
    //myFilezipped->read((char*)&bamtest,bamtestlength);
    bytesread = bgzf_read(bg, bamtest, bamtestlength);
    if(bytesread != bamtestlength){
	cerr<<"Error: GlacParser tried to read "<< bamtestlength <<" bytes but got "<<bytesread<<endl;
	exit(1);
    }
    
    bool startWithBam=true;
    for(unsigned int i=0;i<bamtestlength;i++){
	if(magicNumberBAM[i] != bamtest[i]){
	    startWithBam=false;
	}
    }
    const string magicNumberACF="ACF";
    bool acfFound=true;
    const string magicNumberGLF="GLF";
    bool glfFound=true;

    string toPrintHeader="";
    //unsigned int totalRecords=0;
    string line;    
    istringstream f;
    unsigned char  formattest [5];
    //cerr<<"startWithBam "<<startWithBam<<endl;

    if(!startWithBam){//probably just binary


	bool acftext=true;
	for(unsigned int i=0;i<bamtestlength;i++){
	    firstLineHeadertxtmode+=bamtest[i];
	    if(bamtest[i] != magicNumberACFt[i]){
		acftext=false;
	    }
	}
	if(acftext){

	    stringMode = false;
	    tabixMode  = false;
	    textMode   = true;
	    binMode    = false;

	    glFormat   = false;
	    acFormat   = true;

	}else{
	    
	    bool glftext=true;
	    for(unsigned int i=0;i<bamtestlength;i++){
		if(bamtest[i] != magicNumberGLFt[i]){
		    glftext=false;
		}
	    }
	    if(glftext){
		
		stringMode = false;
		tabixMode  = false;
		textMode   = true;
		binMode    = false;


		glFormat   = true;
		acFormat   = false;				
	    }else{
		cerr<<"Error: GlacParser cannot determine format, first 4 bytes: "<<bamtest[0]<<bamtest[1]<<bamtest[2]<<bamtest[3]<<endl;
		exit(1);		
	    }


	}
	goto findHEADERtext;



	// cerr<<"Error: GlacParser cannot find magic number BAM\1"<<endl;
	// exit(1);

	    
	// formattest [4] = c;
    }else{
	bytesread = bgzf_read(bg, &formattest, 3);
	if(bytesread != 3){
	    cerr<<"Error: GlacParser tried to read "<< 3<<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}

	stringMode=false;
	tabixMode =false;
	textMode  =false;
	binMode   =true;
    
    }


    for(unsigned int i=0;i<3;i++){
	if(formattest[i] != magicNumberACF[i]){
	    acfFound=false;
	}
    }
    
 	
    if(!acfFound){
	//Should find glf
	for(unsigned int i=0;i<3;i++){
	    if(formattest[i] != magicNumberGLF[i]){
		glfFound=false;
	    }
	}

	if(!glfFound){
	    //TODO text mode
	    cerr<<"Error: GlacParser cannot find either acf or glf"<<endl;
	    exit(1);
	}

	
	bytesread = bgzf_read(bg, &sizeBytesGLF, sizeof(sizeBytesGLF));
	if(bytesread != sizeof(sizeBytesGLF)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(sizeBytesGLF) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	if(sizeBytesGLF != 1){
	    cerr<<"Error: GlacParser found  "<< sizeBytesGLF <<" but should 1 "<<endl;
	    exit(1);
	}

	

	char c;
	bytesread = bgzf_read(bg, &c, sizeof(c));
	if(bytesread != sizeof(c)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(c) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}

	if(c!='\1'){
	    cerr<<"Error: GlacParser cannot find magic number GLF\1"<<endl;
	    exit(1);
	}
	glFormat  =true;
	acFormat  =false;
	    	    
    }else{ //ac
	glfFound=false;


	bytesread = bgzf_read(bg, &sizeBytesACF, sizeof(sizeBytesACF));
	if(bytesread != sizeof(sizeBytesACF)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(sizeBytesACF) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	if(sizeBytesACF != 2 && sizeBytesACF != 3){
	    cerr<<"Error: GlacParser found  "<< sizeBytesACF <<" but should 2 or 3 "<<endl;
	    exit(1);
	}



	char c;
	bytesread = bgzf_read(bg, &c, sizeof(c));
	if(bytesread != sizeof(c)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(c) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}

	if(c!='\1'){
	    cerr<<"Error: GlacParser cannot find magic number GLF\1"<<endl;
	    exit(1);
	}
	glFormat  =false;
	acFormat  =true;
	
    }

    //cout<<startWithBam<<"\t"<<acfFound<<"\t"<<glfFound<<"\t"<<int(sizeBytesACF)<<"\t"<<binMode<<endl;    
    uint32_t sizeHeader;
    //myFilezipped.read((char*)&sizeHeader,sizeof(sizeHeader));
    bytesread = bgzf_read(bg, &sizeHeader, sizeof(sizeHeader));
    if(bytesread != sizeof(sizeHeader)){
	cerr<<"Error: GlacParser tried to read "<< sizeof(sizeHeader) <<" bytes but got "<<bytesread<<endl;
	exit(1);
    }
    //cout<<"sh "<<sizeHeader<<endl;
    //exit(1);

    for(uint32_t i=0;i<sizeHeader;i++){
        char toread;
        //myFilezipped.read((char*)&toread,sizeof(char));
	bytesread = bgzf_read(bg, &toread, sizeof(toread));
	if(bytesread != sizeof(toread)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(toread) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	
        toPrintHeader+=toread;
    }


    f.str(toPrintHeader);
    
    while (getline(f, line)) {
        if(strBeginsWith(line,"#SQ")){
            vector<string> tokensf = allTokens(line,'\t');
            chrKnown.push_back(tokensf[1].substr(3));
	    headerSQ+=line+"\n";
        }else{
	    headerNoSQNoDefline+=line+"\n";
	}
	//cout << line << std::endl;
	    // }
	header+=line+"\n";
	if(strBeginsWith(line, "#chr")){
	    defline=line;
	    vector<string> fields=allTokens(line,'\t');
	    if(fields[0] != "#chr")   { cerr<<"Field #1 of header must be #chr ";    exit(1); }
	    if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
	    if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
	    if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
	    if(fields[4] != "anc")    { cerr<<"Field #5 of header must be anc ";     exit(1); }
	    
	    for(unsigned int i=3;i<fields.size();i++){
		populationNames->push_back(fields[i]);
		numberPopulations++;
	    }
	    //header+=line+"\n";
	    
	    break;
	}

    }
    

    //myFilezipped.read((char*)&sizePops,sizeof(sizePops));
    bytesread = bgzf_read(bg, &sizePops, sizeof(sizePops));
    if(bytesread != sizeof(sizePops)){
	cerr<<"Error: GlacParser tried to read "<< sizeof(sizePops) <<" bytes but got "<<bytesread<<endl;
	exit(1);
    }
    
    // if(acfFound)
    // 	sizeRecord = 8+ (2*sizeBytesACF+1)*sizePops;
    // if(glfFound)
    // 	sizeRecord = 8+ (2*sizeBytesGLF+1)*sizePops;
    
    //cout<<sizePops<<endl;
    // exit(1);
    return ;
    
 findHEADERtext:
    //cout<<"header"<<endl;

    //}
    //while(getline ( in,line)){
    //kstring_t * ksstr;
    //cout<<"1"<<endl;
    //while( //bgzf_getline(myFilezipped, '\n', ksstr) != -1){
    line = firstLineHeadertxtmode;
    //line = string(bamtest);
    while(1){
	//line="";
	while(1){
	    char c;
	    bytesread = bgzf_read(bg, &c, sizeof(c));
	    if(bytesread != sizeof(c)){
		cerr<<"Error: GlacParser tried to read "<< sizeof(c) <<" bytes but got "<<bytesread<<endl;
		exit(1);
	    }
	    //cout<<"c"<<c<<"#"<<endl;
	    if(c=='\n'){
		break;
	    }

	    line+=c;
	}
	header+=line+"\n";
	//	cout<<line<<endl;
	// line =string( ks_str(ksstr) );
	// cout<<"line "<<line<<endl;
	//exit(1);

	if(strBeginsWith(line, "#chr")){
	    defline=line;
	    vector<string> fields=allTokens(line,'\t');
	    if(fields[0] != "#chr")   { cerr<<"Field #1 of header must be #chr ";    exit(1); }
	    if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
	    if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
	    if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
	    if(fields[4] != "anc")    { cerr<<"Field #5 of header must be anc ";     exit(1); }

	    for(unsigned int i=3;i<fields.size();i++){
		populationNames->push_back(fields[i]);
		numberPopulations++;
	    }
	    

	    break;
	}
	line="";
    }
    //exit(1);
    // 	//cout<<"line "<<line<<endl;
    // 	if(line[0] == '#'){
    // 	    // cout<<line;
    // 	    if(firstLine){
    // 		if(line != "#MISTAR"){
    // 		    cerr << "Error: GlacParser first line must be #MISTAR found: " << line <<endl;
    // 		    exit(1);	    
    // 		}		
    // 		firstLine=false;
    // 		continue;
    // 	    }

	    
    // 	    if(strBeginsWith(line, "#chr")){
    // 		defline=line;
    // 		vector<string> fields=allTokens(line,'\t');
    // 		if(fields[0] != "#chr")   { cerr<<"Field #1 of header must be #chr ";    exit(1); }
    // 		if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
    // 		if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
    // 		if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
    // 		if(fields[4] != "anc")    { cerr<<"Field #5 of header must be anc ";     exit(1); }

    // 		for(unsigned int i=3;i<fields.size();i++){
    // 		    populationNames->push_back(fields[i]);
    // 		    numberPopulations++;
    // 		}
    // 		header+=line+"\n";

    // 		break;
    // 	    }else{
    // 		header+=line+"\n";
    // 		headerNoDefline+=line+"\n";
    // 	    }
	    
    // 	}else{
    // 	    cerr << "Error: GlacParser cannot get header"  <<endl;
    // 	    exit(1);
    // 	}
    // }
}


string GlacParser::getHeader(string prefix) const{
    vector<string> fields=allTokens(header,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }

    return vectorToString(toreturn,"\n");
}



string GlacParser::getHeaderNoDefline(string prefix) const{
    vector<string> fields=allTokens(headerNoDefline,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }

    return vectorToString(toreturn,"\n");
}

string GlacParser::getHeaderSQ(string prefix) const{
   vector<string> fields=allTokens(headerSQ,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }
    return vectorToString(toreturn,"\n");
}

string GlacParser::getHeaderNoSQNoDefline(string prefix) const{
   vector<string> fields=allTokens(headerNoSQNoDefline,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }
    return vectorToString(toreturn,"\n");
}


string GlacParser::getDefline() const{
    return defline;
}

// void GlacParser::repositionIterator(string chrName,int start,int end){
    
//     if(!tabixMode){
// 	cerr<<"The subroutine repositionIterator can only be called on objects constructed using tabix " <<endl;
// 	exit(1);	
//     }
//     numberOfTimesHasDataWasCalled=-1;

//     rt->repositionIterator(chrName,start,end);
// }

bool GlacParser::hasData(){

    //cout<<"hasData"<<endl;
    if(numberOfTimesHasDataWasCalled!=-1){
	//	cout<<"delete"<<endl;
	//cerr<<"del "<<allRecToReturn<<endl;
	//delete(allRecToReturn->vectorAlleles);
	delete(allRecToReturn);
	numberdel++;
    }else{
	numberOfTimesHasDataWasCalled=0;
    }
    //cout<<"hasDatab1 "<<glFormat<<" "<<binMode<<endl;    
    numberOfTimesHasDataWasCalled++;

    if(tabixMode){

	//AlleleRecords * ar = new AlleleRecords();    
	allRecToReturn                = new AlleleRecords(sizePops,glFormat);
	//cout<<"ar addr"<<allRecToReturn<<endl;
	int result = hts_itr_next(myFilezipped, iter, allRecToReturn, 0);

	    
	//cout<<"result "<<result<<endl;
	if(result ==-1){
	    return false;
	}
	//check if allRecToReturn begin is good
	allRecToReturn->chr = chrKnown[allRecToReturn->chri]; //TODO avoid string copy
	//cout<<"hasData() "<<allRecToReturn->chr<<endl;
	return true;
	//}

    }

    if(binMode){
	numbernew++;
	//cout<<"hasDatab3"<<endl;
	allRecToReturn                = new AlleleRecords(glFormat);
	ssize_t bytesread;


	char ref;
	char alt;
	
	//cout<<"hasDatab2"<<endl;
	//uint16_t 2b
	bytesread = bgzf_read(myFilezipped, &allRecToReturn->chri, sizeof(allRecToReturn->chri));
	if(bytesread == 0){//end of file
	    return false;
	}
	if(bytesread != sizeof(allRecToReturn->chri)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(allRecToReturn->chri) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	//cout<<chrKnown[allRecToReturn->chri]<<endl;
	allRecToReturn->chr = chrKnown[allRecToReturn->chri]; //TODO avoid string copy
	//uint32_t 4b
	bytesread = bgzf_read(myFilezipped, &allRecToReturn->coordinate, sizeof(allRecToReturn->coordinate));
	if(bytesread != sizeof(allRecToReturn->coordinate)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(allRecToReturn->coordinate) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	//cout<<allRecToReturn->coordinate<<endl;
	//char 1b
	bytesread = bgzf_read(myFilezipped, &ref, sizeof(ref));
	if(bytesread != sizeof(ref)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(ref) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	allRecToReturn->ref        =                           "NACGT"[(unsigned char)ref];
	//cout<<"NACGT"[ref]<<endl;
	//char 1b
	bytesread = bgzf_read(myFilezipped, &alt, sizeof(alt));
	if(bytesread != sizeof(alt)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(alt) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	allRecToReturn->alt        =                           "NACGT"[(unsigned char)alt];
	//cout<<"NACGT"[alt]<<endl;
	if(allRecToReturn->ref == allRecToReturn->alt){
	    cerr << "Error: GlacParser the line at "<<allRecToReturn->chr<<":"<<allRecToReturn->coordinate<<" reference allele: "<<allRecToReturn->ref<<" is equal to the alt allele: " << allRecToReturn->alt<<" , exiting"<<endl;
	    exit(1);	   
	}
	
	if(acFormat){
		
	    if(sizeBytesACF ==2 ){ //uint16_t
		
		
		allRecToReturn->vectorAlleles = new vector<SingleAllele>();
		
		//cout<<"sizePops "<<sizePops<<endl;
		
		for(unsigned j=0;j<(sizePops+2);j++){
		    // short refC;
		    // short altC;
		    uint16_t refC;
		    uint16_t altC;

		    char  cpgC;
		    //2b
		    bytesread = bgzf_read(myFilezipped, &refC, sizeof(refC));
		    if(bytesread != sizeof(refC)){
			cerr<<"Error: GlacParser tried to read "<< sizeof(refC) <<" bytes but got "<<bytesread<<endl;
			exit(1);
		    }

		    //2b
		    bytesread = bgzf_read(myFilezipped, &altC, sizeof(altC));
		    if(bytesread != sizeof(altC)){
			cerr<<"Error: GlacParser tried to read "<< sizeof(altC) <<" bytes but got "<<bytesread<<endl;
			exit(1);
		    }
		    //1b
		    bytesread = bgzf_read(myFilezipped, &cpgC, sizeof(cpgC));
		    if(bytesread != sizeof(cpgC)){
			cerr<<"Error: GlacParser tried to read "<< sizeof(cpgC) <<" bytes but got "<<bytesread<<endl;
			exit(1);
		    }

		    //cout<<j<<"\t"<<refC<<"\t"<<altC<<"\t#"<<int(cpgC)<<"#"<<endl;
		    SingleAllele sa (int( refC ),
				     int( altC ),
				     (cpgC == 1) );
		    //cout<<j<<"\t"<<sa<<endl;
		    allRecToReturn->vectorAlleles->push_back(sa);
		    
		}//each pop
		
	    }else{//3 bytes
		exit(1);//to implemenent
	    }
		
	    // cout<<allRecToReturn<<endl;
	    // cout<<*allRecToReturn<<endl;
	    // exit(1);
	    return true;
	}//end if acffound

	// cout<<"sizePops "<<sizePops<<endl;
	// cout<<"gl "<<glFormat<<endl;


	if(glFormat){

	    allRecToReturn->vectorGLs = new vector<SingleGL>();
		
	    // cout<<"sizePops "<<sizePops<<endl;
	    // cout<<glFormat<<endl;
	    //exit(1);
	    
	    for(unsigned j=0;j<(sizePops+2);j++){
		uint8_t rrC;
		uint8_t raC;
		uint8_t aaC;

		uint8_t cpgC;
		    
		bytesread = bgzf_read(myFilezipped, &rrC, sizeof(rrC));
		if(bytesread != sizeof(rrC)){
		    cerr<<"Error: GlacParser tried to read "<< sizeof(rrC) <<" bytes but got "<<bytesread<<endl;
		    exit(1);
		}

		bytesread = bgzf_read(myFilezipped, &raC, sizeof(raC));
		if(bytesread != sizeof(raC)){
		    cerr<<"Error: GlacParser tried to read "<< sizeof(raC) <<" bytes but got "<<bytesread<<endl;
		    exit(1);
		}

		bytesread = bgzf_read(myFilezipped, &aaC, sizeof(aaC));
		if(bytesread != sizeof(aaC)){
		    cerr<<"Error: GlacParser tried to read "<< sizeof(aaC) <<" bytes but got "<<bytesread<<endl;
		    exit(1);
		}
		    
		bytesread = bgzf_read(myFilezipped, &cpgC, sizeof(cpgC));
		if(bytesread != sizeof(cpgC)){
		    cerr<<"Error: GlacParser tried to read "<< sizeof(cpgC) <<" bytes but got "<<bytesread<<endl;
		    exit(1);
		}

		//cout<<j<<"\t"<<refC<<"\t"<<altC<<"\t#"<<int(cpgC)<<"#"<<endl;
		SingleGL gl (rrC,
			     raC,
			     aaC,			     
			     (cpgC == 1) );
		//cout<<j<<"\tgl="<<gl<<"#"<<endl;
		//cout<<"j "<<j<<endl;
		//exit(1);
		allRecToReturn->vectorGLs->push_back(gl);
		
	    }//each pop
	    // cout<<allRecToReturn<<endl;
	    // cout<<*allRecToReturn<<endl;
	    return true;
	}//end glformat
	
    }//end binMode

    if(textMode){
	ssize_t bytesread;
	string line;
	while(1){
	    char c;
	    bytesread = bgzf_read(myFilezipped, &c, sizeof(c));
	    if(bytesread==0 && line.size() == 0){
		return false;
	    }
	    if(bytesread != sizeof(c)){
		cerr<<"Error: GlacParser tried to read "<< sizeof(c) <<" bytes but got "<<bytesread<<endl;
		exit(1);
	    }
	    //cout<<"c"<<c<<"#"<<endl;
	    if(c=='\n'){
		break;
	    }

	    line+=c;
	}
	//cout<<line<<endl;



	allRecToReturn                = new AlleleRecords(glFormat);
	//	cerr<<"new "<<allRecToReturn<<endl;
	//allRecToReturn->vectorAlleles = new vector<SingleAllele>();
	// cout<<"currentline "<<currentline<<endl;
	
	vector<string> fields=allTokens(line,'\t');

	if(fields.size() != (numberPopulations+3)){
	    cerr << "Error: GlacParser the following line should have "<<(numberPopulations+3)<<" fields " << line <<endl;
	    exit(1);	   
	}
	if(fields[2].length() != 3){
	    cerr << "Error: GlacParser the following line " << line <<" does not have 2 comma separated alleles"<<endl;
	    exit(1);	   
	}

	allRecToReturn->chr        =                           fields[0];
	allRecToReturn->coordinate = destringify<unsigned int>(fields[1]);
	allRecToReturn->ref        =                           fields[2][0];
	allRecToReturn->alt        =                           fields[2][2];
	if(allRecToReturn->ref == allRecToReturn->alt){
	    cerr << "Error: GlacParser the following line " << line <<" the reference is equal to the alt allele, exiting"<<endl;
	    exit(1);	   
	}

	if(glFormat){
	    allRecToReturn->vectorGLs = new vector<SingleGL>();

	    for(unsigned int i=3;i<fields.size();i++){
		unsigned int indexComma1=0;
		unsigned int indexComma2=0;		
		unsigned int indexColon=0;
		//cout<<i<<" " <<fields[i]<<endl;
		for(unsigned int k=0;k<fields[i].size();k++){
		    if(fields[i][k]==',' && indexComma1==0)
			indexComma1=k;
		    if(fields[i][k]==',' && indexComma1!=0)
			indexComma2=k;
		    if(fields[i][k]==':')
			indexColon=k;
		}
		
		if(indexComma1 == 0 || indexComma2 == 0 || indexColon == 0 ){
		    cerr << "Error: GlacParser problem with the following line " << line <<" cannot get genotype likelihoods ("<<indexComma1<<","<<indexComma2<<","<<indexColon<<endl;
		    exit(1);	   
		}
		// 		cout<<fields[i].substr(0,indexComma1)<<endl;
		// 		cout<<fields[i].substr(indexComma1+1,indexComma2-indexComma1-1)<<endl;
		// 		cout<<fields[i].substr(indexComma2+1,indexColon-indexComma2-1) <<endl;
		// 		exit(1);
		// // 0,255,0:0
		// //  1   5 7
		// uint8_t t =		destringify<uint8_t>( fields[i].substr(0,indexComma1));
		// cout<<"t "<<t<<endl;
		// cout<<"###"<<endl;
		// cout<<( fields[i].substr(0,indexComma1))<<endl;
		// cout<<( fields[i].substr(indexComma1+1,indexComma2-indexComma1-1))<<endl;
		// cout<<( fields[i].substr(indexComma2+1,indexColon-indexComma2-1))<<endl;
		// uint8_t a= uint8_t( destringify<int>( fields[i].substr(indexComma1+1,indexComma2-indexComma1-1)) );
		// cout<<"i "<< int(a) <<endl;
		// cout<<(fields[i].substr(indexColon+1))<<endl;
	        // cout<<"-------"<<endl;
		
		uint8_t t1 = uint8_t(destringify<int>( fields[i].substr(0,indexComma1)));
		uint8_t t2 = uint8_t(destringify<int>( fields[i].substr(indexComma1+1,indexComma2-indexComma1-1)));
		uint8_t t3 = uint8_t(destringify<int>( fields[i].substr(indexComma2+1,indexColon-indexComma2-1)));
		bool    b  = destringify<bool>(    fields[i].substr(indexColon+1))   ;
		// cout<<indexComma1<<","<<indexComma2<<","<<indexColon<<endl;
		// cout<<t1<<endl;
		// cout<<t2<<endl;
		// cout<<t3<<endl;
		// cout<<b<<endl;
		// cout<<"--------------"<<endl;
		SingleGL gl (t1,
			     t2,
			     t3,
			     b  );
		
		allRecToReturn->vectorGLs->push_back(gl);
	    }
	    
	    if( allRecToReturn->vectorGLs->size() != numberPopulations){
		cerr << "Error: GlacParser problem with the following line " << line <<" number of genotype likelihood read is not "<<numberPopulations<<endl;
		exit(1);	   	    
	    }

	}else{//if acFormat
	    allRecToReturn->vectorAlleles = new vector<SingleAllele>();
	
	    for(unsigned int i=3;i<fields.size();i++){
		unsigned int indexComma=0;
		unsigned int indexColon=0;
		for(unsigned int k=0;k<fields[i].size();k++){
		    if(fields[i][k]==',')
			indexComma=k;
		    if(fields[i][k]==':')
			indexColon=k;
		}
		
		if(indexComma == 0 || indexColon == 0 ){
		    cerr << "Error: GlacParser problem with the following line " << line <<" cannot get allele count"<<endl;
		    exit(1);	   
		}
		
		SingleAllele sa (destringify<int>( fields[i].substr(0,indexComma)),
				 destringify<int>( fields[i].substr(indexComma+1,indexColon)),
			     destringify<bool>(fields[i].substr(indexColon+1))   );
		
		allRecToReturn->vectorAlleles->push_back(sa);
	    }
	    
	    if( allRecToReturn->vectorAlleles->size() != numberPopulations){
		cerr << "Error: GlacParser problem with the following line " << line <<" number of allele count read is not "<<numberPopulations<<endl;
		exit(1);	   	    
	    }
	}//end acFormat
	return true;

	//exit(1);
	//return true;
    }

    //    string line;
    //if(getline ( *myFilezipped,line)){
// if(getNextLine()){
// 	numbernew++;
// 	allRecToReturn                = new AlleleRecords();
// 	//	cerr<<"new "<<allRecToReturn<<endl;
// 	//allRecToReturn->vectorAlleles = new vector<SingleAllele>();
// 	// cout<<"currentline "<<currentline<<endl;
	
// 	vector<string> fields=allTokens(currentline,'\t');

// 	if(fields.size() != (numberPopulations+3)){
// 	    cerr << "Error: GlacParser the following line should have "<<(numberPopulations+3)<<" fields " << currentline <<endl;
// 	    exit(1);	   
// 	}
// 	if(fields[2].length() != 3){
// 	    cerr << "Error: GlacParser the following line " << currentline <<" does not have 2 comma separated alleles"<<endl;
// 	    exit(1);	   
// 	}

// 	allRecToReturn->chr        =                           fields[0];
// 	allRecToReturn->coordinate = destringify<unsigned int>(fields[1]);
// 	allRecToReturn->ref        =                           fields[2][0];
// 	allRecToReturn->alt        =                           fields[2][2];
// 	if(allRecToReturn->ref == allRecToReturn->alt){
// 	    cerr << "Error: GlacParser the following line " << currentline <<" the reference is equal to the alt allele, exiting"<<endl;
// 	    exit(1);	   
// 	}

// 	allRecToReturn->vectorAlleles = new vector<SingleAllele>();
// 	for(unsigned int i=3;i<fields.size();i++){
// 	    unsigned int indexComma=0;
// 	    unsigned int indexColon=0;
// 	    for(unsigned int k=0;k<fields[i].size();k++){
// 		if(fields[i][k]==',')
// 		    indexComma=k;
// 		if(fields[i][k]==':')
// 		    indexColon=k;
// 	    }

// 	    if(indexComma == 0 || indexColon == 0 ){
// 		cerr << "Error: GlacParser problem with the following line " << currentline <<" cannot get allele count"<<endl;
// 		exit(1);	   
// 	    }
	    
// 	    SingleAllele sa (destringify<int>( fields[i].substr(0,indexComma)),
// 			     destringify<int>( fields[i].substr(indexComma+1,indexColon)),
// 			     destringify<bool>(fields[i].substr(indexColon+1))   );

// 	    allRecToReturn->vectorAlleles->push_back(sa);
//  	}

// 	if( allRecToReturn->vectorAlleles->size() != numberPopulations){
// 	    cerr << "Error: GlacParser problem with the following line " << currentline <<" number of allele count read is not "<<numberPopulations<<endl;
// 		exit(1);	   	    
// 	}

// 	return true;

//     }else{//if has no data

// 	if(textMode){
// 	    myFilezipped->close();
// 	}

// 	// else
// 	//     myFile->close();
// 	return false;
//     }
    return false;
}


// bool GlacParser::getNextLine(){
//     if(tabixMode){
// 	return rt->readLine(currentline);
//     }

//     if(textMode){
// 	return getline ( *myFilezipped,currentline);
//     }

//     if(stringMode){
// 	//cout<<dataToReadInd<<"\t"<<dataToRead->size()<<endl;
// 	if(dataToReadInd<dataToRead->size()){
// 	    currentline = dataToRead->at(dataToReadInd++);
// 	    return true;
// 	}else{
// 	    return false;
// 	}
	
//     }
    
//     cerr<<"Invalid state in GlacParser::getNextLine()"<<endl;
//     exit(1);
//     return false;
// }



AlleleRecords * GlacParser::getData(){
    //cout<<"getData"<<allRecToReturn<<endl;
    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;


    return allRecToReturn;
}


bool GlacParser::isACFormat() const{
    return acFormat;
}

bool GlacParser::isGLFormat() const{
    return glFormat;
}

uint32_t GlacParser::getSizePops() const{
    return sizePops;
}

size_t GlacParser::getSizeRecord() const{
    if(acFormat){
	//                    5b record (2b+2b+1)
	size_t sizeSingleAC = (2*sizeBytesACF+1);
	//                  8b base,5b record (2b+2b+1)
	size_t sizeRecord = 8+ (sizeSingleAC)*sizePops;
	return sizeRecord;
    }else{

	if(glFormat){
	    //                     4b record (1+1+1+1)
	    size_t sizeSingleGL = (3*sizeBytesGLF+1);
	    //                  8b base,4b record (1+1+1+1)
	    size_t sizeRecord = 8+ sizeSingleGL*sizePops;
	    return sizeRecord;
	}else{
	    cerr<<"GlacParser: wrong state in getSizeRecord()"<<endl;
	    exit(1);
	}
    }
}

// const vector<string> * GlacParser::getPopulationsNames() const{
//     return populationNames;
// }



