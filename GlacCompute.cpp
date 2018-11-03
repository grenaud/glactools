
#include "GlacCompute.h"


using namespace std;

typedef struct{
    char * buffer;
    unsigned int sizeRecordsRead;
} datachunk;


map<unsigned int, int>       threadID2Rank;
bool isGLF;
char sizeBytesFormat;

//queue< vector<string>  * >  * queueFilesToprocess;
//queue< vector<AlleleRecords>  * >  * queueFilesToprocess;
//queue< char                   * >  * queueFilesToprocess;
queue< datachunk  * >  * queueFilesToprocess;
//queue< char  * >  * queueFilesToprocess;

pthread_mutex_t  mutexTHREADID   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexQueue      = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCounter    = PTHREAD_MUTEX_INITIALIZER;
bool   doneReading;

vector<string> * populationNames;
bool allowUndefined=false;


GlacCompute::GlacCompute(){

}

GlacCompute::~GlacCompute(){

}

template <typename STAT> //type 
void *mainComputationThread(void * argc){

    vector<STAT *> * results =  static_cast<vector<STAT *> *>( argc ) ;

    int   rc;
    // int   stackIndex;
    string freqFileNameToUse;
    int rankThread=0;
    //vector < string >  * dataToUse;
    //vector < AlleleRecords >  * dataToUse;
    datachunk   * dataToUse;

    rc = pthread_mutex_lock(&mutexTHREADID);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
    rankThread = threadID2Rank[*(int *)pthread_self()];
    
    cerr<<"Thread #"<<rankThread <<" is starting"<<endl;

    rc = pthread_mutex_unlock(&mutexTHREADID);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    // cerr<<"Thread #"<<rankThread <<" started and is requesting mutex"<<endl;

    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;

   
    // cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;
    // cerr<<"Thread #"<<(unsigned int)pthread_self() <<" started "<<endl;
    //cout<<"Thread # "<<rankThread<<"TRes "<<results<<endl;

    // cout<<"Thread "<<(unsigned int)pthread_self()<<" taking mutex queue "<<endl;
    if(!queueFilesToprocess->empty()){    
	//cout<<"Thread #"<<rankThread <<" is requesting data"<<endl;
	//cout<<"Thread "<<(unsigned int)pthread_self()<<" taking mutex queue "<<endl;
 	foundData=true;
	dataToUse = queueFilesToprocess->front();
 	queueFilesToprocess->pop();
 	//cerr<<"Thread #"<<rankThread<<" is <<endl;
    }

    //cout<<"Thread #"<<rankThread <<" is requesting data "<<foundData<<" "<<doneReading<<" "<<endl;    
  
    //if(foundData) cout<<"size data "<<dataToUse->size()<<endl;

    if(!foundData){
 	if(doneReading){

	    rc = pthread_mutex_unlock(&mutexQueue);
	    checkResults("pthread_mutex_unlock()\n", rc);

	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
	    return NULL;	
 	}else{
	    //  cout<<"Queue is empty, thread #"<<rankThread<<" will sleep for 5 seconds and wait for data"<<endl;

	    rc = pthread_mutex_unlock(&mutexQueue);
	    checkResults("pthread_mutex_unlock()\n", rc);
	    cerr<<"Queue is empty, thread #"<<rankThread<<" will sleep for 5 seconds, if this happens a lot, consider reducing the number of threads"<<endl;
 	    sleep(2);

 	    goto checkqueue;	   
 	}
    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }


    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////
    //cout<<"Thread #"<<rankThread<<" is starting computations"<<endl;
    // cout<<populationNames<<" "<<vectorToString(*populationNames)<<endl;
    // cout<<"Thread #"<<rankThread<<" is starting computations2"<<endl;
    //unsigned int count=0;
    // for(unsigned int i=0;i<=dataToUse->size();i++){
    // 	//count++;
    // 	cout<<"mcp "<<i<<"\t"<<(dataToUse->at(i))<<endl;
    // }

    //cout<<populationNames<<endl;
    STAT * statComputer = new STAT(populationNames);
    // cout<<"Thread #"<<rankThread <<" addrt stat "<<statComputer<<endl;

    // for(unsigned i=0;i<dataToUse->size();i++){

    // 	//cout<<"Thread #"<<rankThread <<"  "<<i<<endl;

    // 	statComputer->computeStatSingle(&(dataToUse->at(i)),allowUndefined);

    // }
    //cout<<"Thread #"<<rankThread<<" building GP "<<dataToUse->sizeRecordsRead<<" size="<<int(sizeBytesFormat)<<" isGLF "<<isGLF<<endl;

    //cout<<"BUFFERADDR "<<(void*)dataToUse->buffer<<endl;
    GlacParser gp ( dataToUse->buffer,
		    *populationNames,
		    dataToUse->sizeRecordsRead,
		    isGLF,
		    sizeBytesFormat);

    AlleleRecords  * currentRecord;
    unsigned int counterRecords=0;
    //cout<<"Thread #"<<rankThread<<" calling hasData() "<<endl;
    while(gp.hasData()){
	//cout<<"Thread #"<<rankThread<<" before hasData() "<<endl;
    	currentRecord = gp.getData() ;
	//cerr<<"Thread #"<<rankThread<<" "<<*currentRecord<<endl;
	// cout<<"Thread #"<<rankThread<<" after hasData() "<<currentRecord->chri<<":"<<currentRecord->coordinate<<endl;
	if( ((counterRecords%10000)==0) && counterRecords!=0 ){
	    if(rankThread == 2){
		cerr<<"Thread #"<<rankThread<<" is at "<<thousandSeparator(counterRecords)<<endl;
	    }
	}
	counterRecords++;
    	statComputer->computeStatSingle(currentRecord,allowUndefined);

    	//cout<<test->print()<<endl;
    }
    //cout<<"Thread #"<<rankThread<<" DONE hasData() "<<endl;    
    // //TODO set parser
    // GlacParser gp                  ("");    
    // //GlacParser gp                  (dataToUse,*populationNames);    
    // AlleleRecords  * currentRecord;

    // while(gp.hasData()){
    // 	currentRecord = gp.getData() ;

    // 	statComputer->computeStatSingle(currentRecord,allowUndefined);

    // 	//cout<<test->print()<<endl;
    // }
    //cout<<"Thread #"<<rankThread <<" read "<<count<<" records "<<results->size()<<endl;

    //delete(dataToUse);
    //////////////////////////////////////////////////////////////
    //                  END COMPUTATION                         //
    //////////////////////////////////////////////////////////////
    
    delete dataToUse->buffer;
    delete dataToUse;

    //COUNTERS
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);

    //cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;

    results->push_back(statComputer);

    //outputToPrint.push_back(toAddToLog);
    
    //cout<<"Thread #"<<rankThread <<" is re-starting"<<endl;

    rc = pthread_mutex_unlock(&mutexCounter);
    checkResults("pthread_mutex_unlock()\n", rc);


    goto checkqueue;	   


    

    
    //cout<<"Thread "<<rankThread<<" ended "<<endl;
    return NULL;
    
}



template <class STAT> //type 
class parallelP{  

public:
    void launchThreads(const string & filename,int numberOfThreads,int sizeBins,const string & dnaDistMode, const bool performBoot );
};//end class parallelP

//template <class STAT> //type 

// template <typename STAT> 
// vector<STAT> results;

template <class STAT> //type 
void parallelP<STAT>::launchThreads(const string & filename,int numberOfThreads,int sizeBins,const string & dnaDistMode, const bool performBoot ){
    
    doneReading=false;
    //queueFilesToprocess = new queue< vector< string >  * >()  ;
    //queueFilesToprocess = new queue< vector< AlleleRecords >  * >()  ;
    queueFilesToprocess   = new queue< datachunk  * > ();

    //rmd
    // igzstream myFile;
    // myFile.open(filename.c_str(), ios::in);
    // vector<string> lines;

    GlacParser gp (filename);
    vector<string>     chri2chr = gp.getChrKnown();

    if(!gp.isACFormat()){
	cerr<<"GlacCompute: The file "<<filename<<" is not in ACF format"<<endl;
	exit(1);
    }else{
	isGLF=false;
	sizeBytesFormat = int(gp.getSizeOf1DataPoint());
    }
    
    for(unsigned int i=0;i<gp.getPopulationsNames()->size();i++){
	populationNames->push_back(gp.getPopulationsNames()->at(i));
    }

    // vector<string> populationNames;
    //unsigned int numberPopulations;
   
    //allowUndefined  = false;
    //populationNames = mp.getPopulationsNames();

    pthread_mutex_init(&mutexTHREADID,   NULL); 
    pthread_mutex_init(&mutexQueue,      NULL); 
    pthread_mutex_init(&mutexCounter,    NULL);

    pthread_t             thread[numberOfThreads];  
    int                   rc=0;   
    vector<STAT * > * results=new vector<STAT *>();
    //cout<<"res "<<results<<endl;
    //launchThreads(numberOfThreads,*thread);
    for(int i=0;i<numberOfThreads;i++){
        rc = pthread_create(&thread[i], NULL, mainComputationThread<STAT>, results); 
        checkResults("pthread_create()\n", rc); 
    }                                          

    //threads are running here  
    //map<string,int> chr2index;//can be obtained from parser
    //map<string,uint16_t> chri2chr;

    //int chr2indexCurrent=0;
    //uint16_t chriLast=UINT16_MAX;

    //    int indexBin;
    //    int lastBin=-1;
    //int chrBin=-1;

    // AlleleRecords  * currentRecord;
    //vector< string  > * vecForBin;
    // vector< AlleleRecords  > * vecForBin;
    //vector< AlleleRecords  > * vecForBin;
       
    //if (myFile.good()){

    //string line;
    //while ( getline (myFile,line)){
    //AlleleRecords * arr;
    
    //while(gp.hasData()){
    //arr = gp.getData();

    //char * buffer = new char [sizeBlock*gp.getSizeRecord()];
    //char * buffer;
    //int sizeRecordsRead=0;

    datachunk * chunkToAdd      = new datachunk;
    uint32_t coordinate=0;
    uint16_t chri=UINT16_MAX;
    try {
	chunkToAdd->buffer          = new char [sizeBins*gp.getSizeRecord()];
    } catch (std::bad_alloc&) {
	cerr<<"Could not allocate a buffer of size "<<(sizeBins*gp.getSizeRecord())<<endl;
	exit(1);
    }

    //cout<<"CREATING BUFFER1 "<<(void*)chunkToAdd->buffer<<endl;
    
    chunkToAdd->sizeRecordsRead = 0;
    // cout<<"chunkToAdd1 "<<chunkToAdd<<endl;
    // cout<<"chunkToAdd buffer #"<<(void *) chunkToAdd->buffer<<"#"<<endl;
    // cout<<"buffer size "<<(sizeBins*gp.getSizeRecord())<<endl;
    // cerr<<"reading block "<<sizeBins<<endl;
    // cout<<"chunkToAdd2 "<<chunkToAdd<<" "<<&chunkToAdd<<endl;

    
    bool rbdReturn = gp.readBlockData(chunkToAdd->buffer,sizeBins,&chunkToAdd->sizeRecordsRead,&chri,&coordinate);
    // cout<<"rbdReturn "<<rbdReturn<<endl;
    // cout<<"chunkToAdd3 "<<chunkToAdd<<" "<<&chunkToAdd<<endl;

    while( rbdReturn ){
	//cout<<*arr<<endl;
	cerr<<"GlacCompute reading new bin, currently  "<<chri2chr[chri]<<":"<<thousandSeparator(coordinate)<<" size="<<chunkToAdd->sizeRecordsRead<<endl;
	
	int rc = pthread_mutex_lock(&mutexQueue);
	checkResults("pthread_mutex_lock()\n", rc);
	
		
	bool needToAskMutex=false;
	//add old
	while( int(queueFilesToprocess->size()) > numberOfThreads ){
	    needToAskMutex=true;
	    //unlock mutex
	    rc = pthread_mutex_unlock(&mutexQueue);
	    checkResults("pthread_mutex_unlock()\n", rc);
	    
	    cerr<<"queue is full, sleeping for 4 seconds"<<endl;
	    sleep(4);
	}

	if(needToAskMutex){
	    rc = pthread_mutex_lock(&mutexQueue);
	    checkResults("pthread_mutex_lock()\n", rc);
	}
		
	queueFilesToprocess->push(chunkToAdd);

	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
	
	chunkToAdd                  = new datachunk;
	chunkToAdd->buffer          = new char [sizeBins*gp.getSizeRecord()];
	//cout<<"CREATING BUFFER2 "<<(void*)chunkToAdd->buffer<<endl;
	chunkToAdd->sizeRecordsRead = 0;
	

	//cout<<"p1#"<<arr->chr<<"\t"<<arr->coordinate<<endl;
	//vecForBin->push_back(*arr);
	//vecForBin->push_back(*arr);

	//cout<<"p2#"<<arr->chr<<"\t"<<arr->coordinate<<endl;
	//vecForBin->push_back(line);
	rbdReturn = gp.readBlockData(chunkToAdd->buffer,sizeBins,&chunkToAdd->sizeRecordsRead,&chri,&coordinate);
    }//done reading

    cerr<<"done reading, adding final chunk at "<<chri2chr[chri]<<":"<<thousandSeparator(coordinate)<<" size="<<chunkToAdd->sizeRecordsRead<<endl;
    // cout<<"chunkToAdd "<<chunkToAdd<<endl;
    // cout<<"chri "<<chri<<endl;
    // cout<<"coordinate "<<coordinate<<endl;
    // cout<<"done reading  "<< chunkToAdd->sizeRecordsRead<<endl;


    //adding the last chunk
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);
	
		
    bool needToAskMutex=false;
    //add old
    while( int(queueFilesToprocess->size()) > numberOfThreads ){
	needToAskMutex=true;
	//unlock mutex
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
	    
	cerr<<"queue is full, sleeping for 4 seconds"<<endl;
	sleep(4);
    }

    if(needToAskMutex){
	rc = pthread_mutex_lock(&mutexQueue);
	checkResults("pthread_mutex_lock()\n", rc);
    }
    
    if(chunkToAdd->sizeRecordsRead>0)
	queueFilesToprocess->push(chunkToAdd);

    rc = pthread_mutex_unlock(&mutexQueue);
    checkResults("pthread_mutex_unlock()\n", rc);
	

    //queueFilesToprocess->push(vecForBin);//adding last bin
    //queueFilesToprocess->push(buffer);//adding last bin
    //queueFilesToprocess->push(chunkToAdd);//adding last bin

    //myFile.close();

    // }else{
    // 	cerr << "Unable to open file "<<filename<<endl;
    // 	exit(1);
    // }
    


    doneReading=true;
    //waiting for threads to finish
    for (int i=0; i <numberOfThreads; ++i) {
        rc = pthread_join(thread[i], NULL); 
        checkResults("pthread_join()\n", rc);
    }
    //cout<<"ALL DONE1"<<endl;
    pthread_mutex_destroy(&mutexQueue);   
    pthread_mutex_destroy(&mutexCounter);
    pthread_mutex_destroy(&mutexTHREADID);
    //cout<<"ALL DONE2"<<endl;



    //DO jacknifing
    if(results->size()>1){
	//cout << "VEC "<<vectorToString( *((*results->at(0)).populationNames) )<<endl;
	STAT * allResults=new STAT((*results->at(0)));
	// cout << "VEC "<<vectorToString( *((*results->at(0)).populationNames) )<<endl;
	// cout<<"done all1"<<endl;
	// cout<<allResults->print()<<endl;
	// cout<<"done all2"<<endl;
	vector<STAT * > * jacknife=new vector<STAT *>();

	for (unsigned int i=1; i<results->size() ; ++i) {    
	    //cout<<"add\t"<<i<<endl;	
	    (*allResults)+=(*results->at(i));
	    // //cout<<"done dd1"<<endl;
	    // cout<<allResults->print()<<endl;
	    // cout<<"done dd2"<<endl;
	}

	for (unsigned int i=0; i<results->size() ; ++i) {    
	    // cout<<i<<"\n#####\n"<<endl;
	    //	cout<<results->at(i)<<endl;
	    cout<<"---------------------------"<<endl;
	    cout<<results->at(i)->print()<<endl;
	    //	cout<<i<<"\n#####\n"<<endl;
	}
	
	if( performBoot ){
	    for (unsigned int i=0; i<results->size() ; ++i) {    
		STAT * test =new STAT (*allResults); //creating a copy
		*test-=(*results->at(i)); //removing ith block
		jacknife->push_back(test); 
		//cout<<"ji "<<i<<endl<<test->print()<<endl;
	    }
	    
	    cout<<"---------------ALL---------------"<<endl;
	    //cerr<<"done all1"<<endl;
	    //COMMENT allResults contains a matrix of AvgCoaResults
	    //add a method for jacknife in allResults
	    cout<<allResults->printWithBootstraps(jacknife,dnaDistMode);
	    //cout<<allResults->print()<<endl;
	    //cout<<"done all2"<<endl;
	}
    }else{
	if( !performBoot ){
	    cout<<"---------------------------"<<endl;
	    cout<<results->at(0)->print()<<endl;
	}else{//need to perform boot
	    cerr<<"GlacCompute: There is only a single bin, cannot perform a jacknife, please run with more data"<<endl;
	}
    }
    pthread_exit(NULL); 
    //cout<<"ALL DONE3"<<endl;
}

string GlacCompute::usage() const{
    string usage=string("glactools")+" compute [options]  <ACF file>"+
                                    "\nThis program computers summary stats on ACF files\n\n"+
	
	                            "Options:\n"+ 
	"\t"+"-p [stats]"  +"\t\t" +"Statistics to use:\n"+
         			      "\t"+"    paircoacompute"+"\tTo compute pairwise average coalescence\n"+
	//			      "\t"+"    fst"+"\t\t\tTo compute pairwise Fst (Weir and Cockerham's 1984)\n"+
			    "\t"+"    dstat"+"\t\tTo compute triple-wise D-statistics\n"+
        	            "\t"+"    dist"+"\t\tTo compute simple pairwise distance\n"+
	                    "\t"+"    f3  "+"\t\tF3 stats (Pickrell implementation)\n"+
	                    "\t"+"    f2  "+"\t\tF2 stats (Pickrell implementation)\n"+

	"\n"+
                       	//"\t"+"    f2"+"\t\tF2 stats (Pickrell implementation)\n\n"+
			    //   "\t"+"    For the \"dist\" mode, specify the model:\n"+
			    //   "\t"+"    "+"\t--model [model]\tUse this model for DNA distance\n"+
                            //   "\n"+
	                    //   "\t"+"    "+"\tnone\t"+           "all mutations with equal footing (default)\n"+
			    //   "\t"+"    "+"\tJC69\t"+		"Jukes Cantor 1969\n"+
			    //   "\t"+"    "+"\tK80\t"+		"Kimura 1980\n"+
                            // //"\t"+"    "+"\t\t\tHKY85\t"+		"Hasegawa, Kishino and Yano 1985\n"+
                              "\n"+

	"\t"+"-u"  +"\t\t\t" +"Allow undefined sites (0,0) for certain individuals/pops, can cause ascertainment bias (Default: "+stringify(allowUndefined)+")\n"+
	"\t"+"-t [threads]"  +"\t\t" +"Threads to use (Default: "+stringify(numberOfThreads)+")\n"+
        "\t"+"-s [size bin]" +"\t\t" +"Size of bins (Default: "+stringify(sizeBins)+")\n"+
        "\t"+"--noboot"      +"\t\t" +"Skip the bootstrap, recommended for large datasets (Default: "+booleanAsString(!performBoot)+")\n"+
        "\t"+"--boot"        +"\t\t" +"If you have multiple chunks computed using --noboot, use this option to combine them\n"+
        "\t"+"     "        +"\t\t"  +"and just perform the bootstaps, you need to specify the option files as arguments\n"+
        "\t"+"     "        +"\t\t"  +"(Default: "+booleanAsString(justBoot)+")\n"+
	"\n";


    return usage;
}



template <class STAT> //type 
void GlacCompute::bootFromResults(vector<string> * arguments,STAT * stattouse ){
    ////////////////////////////////
    //READING PREVIOUS RESULTS   ///
    ////////////////////////////////
    vector<STAT * > * results=new vector<STAT *>();
    //cerr<<"bootFromResults"<<endl;
    for(unsigned int i=0;i<arguments->size();i++){
	igzstream myfileResults;
	STAT * statComputer = new STAT();
	myfileResults.open(arguments->at(i).c_str(), ios::in);
	if (myfileResults.good()){
	    //fix the operators
	    // while (!myfileResults.eof()){
	    // 	myfileResults >> statComputer ;

	    // }

	    cerr<<"Reading file "<<arguments->at(i)<<endl;
	    string strResults;
	    string lineResults;
	    while( getline(myfileResults,lineResults)){
	    	strResults+=lineResults+"\n";
	    }
	    //istringstream in (strResults);
	    statComputer->read(strResults);
	    results->push_back(statComputer);
	    //cout<<*statComputer<<endl;
	}else{
	    cerr<<"Cannot open file "<<arguments->at(i)<<endl;
	    exit(1);
	}
    }
    //    exit(1);

    cerr<<"reading done, starting jacknifing"<<endl;
    //DO jacknifing
    if(results->size()>1){
	//cout << "VEC "<<vectorToString( *((*results->at(0)).populationNames) )<<endl;
	STAT * allResults=new STAT((*results->at(0)));
	// cout << "VEC "<<vectorToString( *((*results->at(0)).populationNames) )<<endl;
	// cout<<"done all1"<<endl;
	// cout<<allResults->print()<<endl;
	// cout<<"done all2"<<endl;
	vector<STAT * > * jacknife=new vector<STAT *>();

	for (unsigned int i=1; i<results->size() ; ++i) {    
	    //cout<<"add\t"<<i<<endl;	
	    (*allResults)+=(*results->at(i));
	    // //cout<<"done dd1"<<endl;
	    // cout<<allResults->print()<<endl;
	    // cout<<"done dd2"<<endl;
	}

	// for (unsigned int i=0; i<results->size() ; ++i) {    
	//     // cout<<i<<"\n#####\n"<<endl;
	//     //	cout<<results->at(i)<<endl;
	//     cout<<"---------------------------"<<endl;
	//     cout<<results->at(i)->print()<<endl;
	//     //	cout<<i<<"\n#####\n"<<endl;
	// }
	
	if( performBoot ){
	    for (unsigned int i=0; i<results->size() ; ++i) {    
		cerr<<"jackknife #"<<i<<" of "<<results->size()<<endl;
		STAT * test =new STAT (*allResults); //creating a copy
		*test-=(*results->at(i)); //removing ith block
		jacknife->push_back(test); 
		//cout<<"ji "<<i<<endl<<test->print()<<endl;
	    }
	    
	    cout<<"---------------ALL---------------"<<endl;
	    //cerr<<"done all1"<<endl;
	    //COMMENT allResults contains a matrix of AvgCoaResults
	    //add a method for jacknife in allResults
	    cout<<allResults->printWithBootstraps(jacknife,dnaDistMode);
	    //cout<<allResults->print()<<endl;
	    //cerr<<"done all2"<<endl;
	}

	delete(jacknife);
	delete(allResults);
    }else{
	if( !performBoot ){
	    cout<<"---------------------------"<<endl;
	    cout<<results->at(0)->print()<<endl;
	}else{//need to perform boot
	    cerr<<"GlacCompute: There is only a single bin, cannot perform a jacknife, please run with more data"<<endl;
	}
    }
    cerr<<"jacknifing done"<<endl;
    delete(results);
}

int GlacCompute::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    populationNames = new 	vector<string>  ();
    //cout<<"run()"<<endl;
    //cout<<populationNames<<" "<<vectorToString(*populationNames)<<endl;

    bool dnaDistModeSpecified=false;
    dnaDistMode              ="none";
    int lastOpt=1;

    for(int i=1;i<(argc-1);i++){ 
	//cerr<<i<<" "<<string(argv[i])<<endl;
	if((string(argv[i]) == "-")  ){
            lastOpt=i;
            break;          
        }
	
	if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }

	
        if(string(argv[i]) == "-t" ){
            numberOfThreads=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-u" ){
            allowUndefined=true;
            continue;
        }

        if(string(argv[i]) == "--noboot" ){
	    performBoot=false;
            continue;
        }

        if(string(argv[i]) == "--boot" ){
	    justBoot=true;
            continue;
        }

        if(string(argv[i]) == "-s" ){
	    sizeBins=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-p" ){
	    program=string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-p" ){
	    program=string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "--model" ){
	    dnaDistMode=string(argv[i+1]);
	    dnaDistModeSpecified=true;
            i++;
            continue;
        }

        cerr<<"Wrong option "<<argv[i]<<endl;
        return 1;

    }
    
    if( dnaDistModeSpecified ){
	if(program != "dist"){
	    cerr<<"GlacCompute: ERROR: Cannot specify --model if the option \"dist\" is not used "<<endl;
	    return 1;	    
	}
    }

    if(justBoot){//just perform bootstaps
	if(program == "dstat"){
	    //parallelP<SumStatD> pToRun;
	    //pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);
	    
	    vector<string> * arguments=new vector<string>();
	    for(int i=lastOpt;i<argc;i++){
		//cerr<<string(argv[i])<<endl;
		arguments->push_back( string(argv[i]) );
	    }
	    
	    SumStatD * st=new SumStatD();

	    bootFromResults(arguments,st);

	    delete(st);
	    delete(arguments);
	    
	    // for(int i=1;i<(argc-1);i++){ 
	    // 	if((string(argv[i]) == "-")  ){
	    // 	    lastOpt=i;
	    // 	    break;          
	    // 	}
		
	}else{
	    cerr<<"GlacCompute: to implement (coming soon) "<<endl;
	    return 1;	    

	}
    }else{
	if(program == "paircoacompute"){
	    parallelP<SumStatAvgCoa> pToRun;
	    pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);
	}else{
	    if(program == "dstat"){
		parallelP<SumStatD> pToRun;
		pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);
	    }else{
		if(program == "fst"){
		    parallelP<SumStatFst> pToRun;
		    pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);
		}else{
		    if(program == "dist"){
			parallelP<SumStatDist> pToRun;
			pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);	    
		    }else{
			if(program == "f3"){
			    parallelP<SumStatF3> pToRun;
			    pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);	    
			}else{
			    if(program == "f2"){
				parallelP<SumStatF2> pToRun;
				pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,dnaDistMode,performBoot);	    
			    }else{
				cerr<<"Wrong program "<<program<<endl;
				return 1;
			    }
			}
		    }
		}
	    }
	}
    }	
    return 0;
}



// #ifdef OLD
// 	exit(1);
// 	if( chri !=  chriLast ){  //=UINT16_MAX)
// 	    //chr2index[ chrS ]  = (chr2indexCurrent++);//adding
// 	    chriLast = chri;
// 	    // if(lastBin != -1)
// 	    // 	currentBin = lastBin+1+currentBin;
// 	    if(chrBin == -1){
// 		chrBin = 0;
// 	    }else{
// 		chrBin = lastBin +1; //a step above the last bin
// 	    }
	    
// 	    cerr<<"new chr found, processing chr #: "<<chri2chr[ chri ]<<endl;
// 	    //(chrS)<<endl;
// 	}
	
// 	int currentBin = chrBin+(coordinate/sizeBins);	
// 	cout<<"GlacCompute chri "<<lastBin<<" "<<currentBin<<endl;

// 	if(lastBin != currentBin){
// 	    //cout<<"2: "<<arr->chr<<"\tc="<<arr->coordinate<<"\tbin="<<(currentBin)<<endl;

// 	    //cout<<"new bin"<<endl;
// 	    if( (currentBin%10)==0){
// 		//cerr<<"processing : "<<chrBin<<":"<<coordSUI<<endl;
// 		cerr<<"processing chr #: "<<(chri2chr[ chri ])<<":"<<coordinate<<endl;
// 	    }

// 	    if(lastBin == -1){ //first bin
// 		cout<<"new bin first"<<endl;
// 		//re-init
// 		//vecForBin =  new vector< AlleleRecords > ();				
// 		//vecForBin->reserve(sizeBins);
// 		//buffer = new char [sizeBlock*gp.getSizeRecord()];
// 		datachunk * chunkToAdd      = new datachunk;
// 		chunkToAdd->buffer          = new char [sizeBins*gp.getSizeRecord()];
// 		cout<<"CREATING BUFFER2 "<<(void*)chunkToAdd->buffer<<endl;
// 		chunkToAdd->sizeRecordsRead = 0;
		
// 		//vecForBin =  new vector< string > ();		

// 		lastBin=currentBin;
// 	    }else{

// 		int rc = pthread_mutex_lock(&mutexQueue);
// 		checkResults("pthread_mutex_lock()\n", rc);
// 		//cout<<"new bin old\t"<<int(queueFilesToprocess->size())<<"\tresadr\t"<<results<<"\tressize\t"<<results->size()<<endl;
		
// 		bool needToAskMutex=false;
// 		//add old
// 		while( int(queueFilesToprocess->size()) > numberOfThreads ){
// 		    needToAskMutex=true;
// 		    //unlock mutex
// 		    rc = pthread_mutex_unlock(&mutexQueue);
// 		    checkResults("pthread_mutex_unlock()\n", rc);

// 		    //cout<<"Queue is full main threads will sleep for 10 seconds and wait for threads to finish"<<endl;
// 		    sleep(4);
// 		}

// 		if(needToAskMutex){
// 		    rc = pthread_mutex_lock(&mutexQueue);
// 		    checkResults("pthread_mutex_lock()\n", rc);
// 		}
		
// 		//queueFilesToprocess->push(vecForBin);
// 		//buffer = new char [sizeBlock*gp.getSizeRecord()];
// 		queueFilesToprocess->push(chunkToAdd);

// 		rc = pthread_mutex_unlock(&mutexQueue);
// 		checkResults("pthread_mutex_unlock()\n", rc);
		
// 		//cout<<"pushing "<<vecForBin<<endl;
// 		// for(unsigned int j=0;j<vecForBin->size();j++){
// 		//     cout<<"v "<<(vecForBin->at(j))<<endl;
// 		// }
// 		//re-init
// 		//vecForBin =  new vector< string  > ();		
// 		//vecForBin =  new vector< AlleleRecords  > ();
// 		//vecForBin->reserve(sizeBins);
// 		//buffer = new char [sizeBlock*gp.getSizeRecord()];
// 		datachunk * chunkToAdd      = new datachunk;
// 		chunkToAdd->buffer          = new char [sizeBins*gp.getSizeRecord()];
// 		cout<<"CREATING BUFFER3 "<<(void*)chunkToAdd->buffer<<endl;
// 		chunkToAdd->sizeRecordsRead = 0;


// 		lastBin=currentBin;
// 		//
// 	    }
// 	}else{
// 	    //cout<<"old bin"<<endl;
// 	}
// #endif	
