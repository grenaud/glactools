/*
 * GlacParser
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GlacParser.h"



// GlacParser::GlacParser(string file,string indexForFile,string chrName,int start,int end){
//     //for rand() ins SingleAllele
//     timeval time;
//     gettimeofday(&time, NULL);
//     srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );


//     rt =new ReadTabix (file,indexForFile,chrName,start,end);
//     string headertemp=rt->getHeader();

//     istringstream is(headertemp);
  
//     //reading the header
//     numberPopulations=0;
//     populationNames=new vector<string>();

//     parseHeader(is);

//     if(defline.empty()){
// 	cerr << "Error: GlacParser cannot get definition line"  <<endl;
// 	exit(1);
//     }

//     numberOfTimesHasDataWasCalled=-1;
//     tabixMode  = true;
//     textMode   = false;
//     stringMode = false;
// }


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

    bool isbgzip=bgzf_is_bgzf(filename.c_str());

    myFilezipped = bgzf_open(filename.c_str(), "r");

    if(myFilezipped == NULL){
        cerr<<"Error: GlacParser failed to open file "<< filename <<endl;
        exit (1);
    }

    if(isbgzip){
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
    //cerr<<"destructor GlacParser"<<endl;
    if(numberOfTimesHasDataWasCalled == 0){//last called was getData
	//delete(allRecToReturn->vectorAlleles);
	delete(allRecToReturn);
	numberdel++;
    }

    if(tabixMode){
	delete rt; //calling the destructor
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
    unsigned int totalRecords=0;
    string line;    
    istringstream f;
    unsigned char  formattest [5];
    if(!startWithBam){//probably just binary
	//cout<<startWithBam<<endl;


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
        }//else{
	cout << line << std::endl;
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
	    header+=line+"\n";
	    
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


string GlacParser::getHeader(string prefix){
    vector<string> fields=allTokens(header,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }

    return vectorToString(toreturn,"\n");
}



string GlacParser::getHeaderNoDefline(string prefix){
    vector<string> fields=allTokens(headerNoDefline,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }

    return vectorToString(toreturn,"\n");
}




string GlacParser::getDefline(){
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

    if(binMode){
	numbernew++;
	//cout<<"hasDatab3"<<endl;
	allRecToReturn                = new AlleleRecords(glFormat);
	ssize_t bytesread;


	char ref;
	char alt;
	
	//cout<<"hasDatab2"<<endl;
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
	bytesread = bgzf_read(myFilezipped, &allRecToReturn->coordinate, sizeof(allRecToReturn->coordinate));
	if(bytesread != sizeof(allRecToReturn->coordinate)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(allRecToReturn->coordinate) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	//cout<<allRecToReturn->coordinate<<endl;

	bytesread = bgzf_read(myFilezipped, &ref, sizeof(ref));
	if(bytesread != sizeof(ref)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(ref) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	allRecToReturn->ref        =                           "NACGT"[ref];
	//cout<<"NACGT"[ref]<<endl;
	bytesread = bgzf_read(myFilezipped, &alt, sizeof(alt));
	if(bytesread != sizeof(alt)){
	    cerr<<"Error: GlacParser tried to read "<< sizeof(alt) <<" bytes but got "<<bytesread<<endl;
	    exit(1);
	}
	allRecToReturn->alt        =                           "NACGT"[alt];
	//cout<<"NACGT"[alt]<<endl;
	if(allRecToReturn->ref == allRecToReturn->alt){
	    cerr << "Error: GlacParser the following line " << currentline <<" the reference is equal to the alt allele, exiting"<<endl;
	    exit(1);	   
	}
	
	if(acFormat){
		
	    if(sizeBytesACF ==2 ){
		
		
		allRecToReturn->vectorAlleles = new vector<SingleAllele>();
		
		//cout<<"sizePops "<<sizePops<<endl;
		
		for(unsigned j=0;j<(sizePops+2);j++){
		    short refC;
		    short altC;
		    char  cpgC;
		    
		    bytesread = bgzf_read(myFilezipped, &refC, sizeof(refC));
		    if(bytesread != sizeof(refC)){
			cerr<<"Error: GlacParser tried to read "<< sizeof(refC) <<" bytes but got "<<bytesread<<endl;
			exit(1);
		    }

		    bytesread = bgzf_read(myFilezipped, &altC, sizeof(altC));
		    if(bytesread != sizeof(altC)){
			cerr<<"Error: GlacParser tried to read "<< sizeof(altC) <<" bytes but got "<<bytesread<<endl;
			exit(1);
		    }
		    
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
	cout<<line<<endl;



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
		for(unsigned int k=0;k<fields[i].size();k++){
		    if(fields[i][k]==',' && indexComma1==0)
			indexComma1=k;
		    if(fields[i][k]==',' && indexComma1!=0)
			indexComma1=k;
		    if(fields[i][k]==':')
			indexColon=k;
		}
		
		if(indexComma1 == 0 || indexComma2 == 0 || indexColon == 0 ){
		    cerr << "Error: GlacParser problem with the following line " << line <<" cannot get genotype likelihoods"<<endl;
		    exit(1);	   
		}
		
		SingleGL gl (destringify<uint8_t>( fields[i].substr(0,indexComma1)),
			     destringify<uint8_t>( fields[i].substr(indexComma1+1,indexComma2)),
			     destringify<uint8_t>( fields[i].substr(indexComma2+1,indexColon)),
			     destringify<bool>(fields[i].substr(indexColon+1))   );
		
		allRecToReturn->vectorGLs->push_back(gl);
	    }
	    
	    if( allRecToReturn->vectorGLs->size() != numberPopulations){
		cerr << "Error: GlacParser problem with the following line " << line <<" number of genotype likelihood read is not "<<numberPopulations<<endl;
		exit(1);	   	    
	    }

	}else{
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
	}
	return true;

	exit(1);
	return true;
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





// const vector<string> * GlacParser::getPopulationsNames() const{
//     return populationNames;
// }



