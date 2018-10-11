
#include "GlacWindows.h"


using namespace std;


GlacWindows::GlacWindows(){

}

GlacWindows::~GlacWindows(){

}

string GlacWindows::usage() const{
    string usage=string("glactools")+" windows  <ACF or GLF file>"+
	"\nThis program will print set of windows along the genome\n"+
	"Options:\n"+
	//"\t\t"+"--allowsexchr" +"\t\t\t\t"+"Allow sites to be generated on X and Y (default: "+boolStringify(allowSexChr)+")\n";
	"\n\tLoci size option:\n"+
	"\t\t"+"--chunk [size]"     +"\t\t\t\t"+"Size of contiguous genomic region to take (default: none)\n"+
	"\t\t"+"--size  [file]"     +"\t\t\t\t"+"Instead of specifying --chunk for --random, read this file with the sizes\n"+
	"\t\t                 "     +"\t\t\t"+"(one per line, no header, just integers)\n"+
	"\t\t"+"--overlap [size]"   +"\t\t\t"+  "Size of overlap between windows (for overlap mode) (default: "+stringify(overlapBetweenWindows)+")\n"+
  	"\n\tLoci selection mode, select either:\n"+
	// "\t\t"+"--region [file]"    +"\t\t\t\t"   +"Read [file] in BED format to select which loci to use\n"+
	"\t\t"+"--random [number]"  +"\t\t\t"     +"Select [number] random loci from the genome\n"+
	"\t\t"+"--loci"             +"\t\t\t\t\t" +"Select multiple loci from the genome\n"+
	"\t\t"+"--locichr [chr]"    +"\t\t\t\t"   +"Select multiple loci from a chromosome\n"+
	"\t\t"+"--chr [chr]"        +"\t\t\t\t"   +"Run on an entire chromosome\n"+
	"\t\t"+"--wide"             +"\t\t\t\t\t" +"Run on the entire genome\n"+
      "\n\tOutput options:\n"+
      "\t\t"+"--bed"  +"\t\t\t"     +"Produce output in bed format (default: "+boolStringify(inBedFormat)+")\n\n";
    return usage;
}

int GlacWindows::run(int argc, char *argv[]){

    //int main (int argc, char *argv[]) {


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }


    int lastOpt=1;

    for(int i=1;i<(argc-1);i++){ 
        //cout<<i<<"\t"<<string(argv[i])<<endl;
        if((string(argv[i]) == "-")  ){
            lastOpt=i;
            break;          
        }

        if(string(argv[i]) == "--bed" ){
	    inBedFormat=true;
            continue;
        }

	//loci size
        if(string(argv[i]) == "--chunk" ){
	    bpToExtract=destringify<int>(argv[i+1]);
	    i++;
	    specifiedChunk=true;
            continue;
        }

        if(string(argv[i]) == "--size" ){
	    sizeFile=string(argv[i+1]);
	    specifiedSizeFile=true;
	    i++;
            continue;
        }

	if(string(argv[i]) == "--overlap" ){
	    overlapBetweenWindows=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }


	if(string(argv[i]) == "--allowsexchr" ){
	    allowSexChr=true;
	    continue;
	}

	//BEGIN loci selection
        if(string(argv[i]) == "--random" ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    amountOfGoodSitesTARGET=destringify<unsigned int>(argv[i+1]);
	    lociSelection='r';
	    i++;
            continue;
        }

        if(string(argv[i]) == "--loci" ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }

	    lociSelection='o';
            continue;
        }

        if(string(argv[i]) == "--locichr" ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    chrToUse=string(argv[i+1]);
	    i++;	    
	    lociSelection='c';
            continue;
        }

        // if(string(argv[i]) == "--region" ){
	//     if(lociSelection != 'x'){
	// 	cerr << "Choose only one loci selection option"<<endl;
	// 	return 1;       
	//     }
	//     bedFileWithRegions=string(argv[i+1]);
	//     i++;	    
	//     lociSelection='f';
        //     continue;
        // }


        if(string(argv[i]) == "--chr" ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    chrToUse=string(argv[i+1]);
	    i++;	    
	    lociSelection='v'; 
            continue;
        }

	if(string(argv[i]) == "--wide" ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    lociSelection='w';
            continue;
        }


        cerr<<"Error unknown option "<<argv[i]<<endl;
        return 1;

    }




    if(specifiedChunk && specifiedSizeFile){
	cerr << "Error, cannot specify both --size and --chunk"<<endl;
	return 1;       
    }

    if(lociSelection != 'r' &&
       specifiedSizeFile ){ //random
	cerr << "The option --size is only for random genomic region"<<endl;
	return 1;       	
    }

    if(lociSelection != 'v' &&  //entire chr
       lociSelection != 'w' ){	 //genome wide
	if(!specifiedChunk && !(lociSelection=='r' && specifiedSizeFile) ){
	    cerr << "Need to specify the size (with --chunk) of the size of contiguous genomic region to take"<<endl;
	    return 1;       
	}

    }


    vector<int> sizeOfBp;

    if(specifiedSizeFile){

	string line;
	igzstream myFile;
	string filename = string(sizeFile);
	myFile.open(filename.c_str(), ios::in);
	
	if (myFile.good()){
	    while ( getline (myFile,line)){
		int toadd=destringify<int>(line);
		if(toadd<0){
		    cerr << "Cannot have negative values in file for --size line="<<line<<endl;
		    return 1;     
		}
		    
		sizeOfBp.push_back(toadd);
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<filename<<endl;
	    return 1;
	}

	if(sizeOfBp.empty()){
	    cerr << "Error:  file "<<filename<<" does not contain any size"<<endl;
	    return 1;
	}
	
    }


    string glacfile  = string(argv[argc-1]);


    
    GlacParser gp (glacfile);
    string sqLines = gp.getHeaderSQ();

    // vector<string> sq = allTokens(linesSQ,'\n');
    // for(unsigned int i=0;i<sq.size();i++){
    // 	vector<string> sql = allTokens(sq[i],'\t');
    // 	cout<<sql[2]<<endl;	
    // }
	

    RandomGenomicCoord rgc (sqLines,allowSexChr);
    GenomicWindows     rw  (sqLines,allowSexChr);
    vector<GenomicRange> v;
    unsigned int indexOfLoci=0;
    switch(lociSelection){
    case 'r': //random loci
	
	while(indexOfLoci<amountOfGoodSitesTARGET){
	    if(specifiedSizeFile){
		//srand is called in RandomGenomicCoord
		unsigned int toReturn = (unsigned int)(mrand48());
		bpToExtract=sizeOfBp[toReturn%sizeOfBp.size()];
	    }
	    // v.push_back(rgc.getRandomGenomic(bpToExtract));
	    indexOfLoci++;
	    if(!inBedFormat)
		cout<<rgc.getRandomGenomic(bpToExtract)<<endl;
	    else
		cout<<rgc.getRandomGenomic(bpToExtract).asBed()<<endl;
	}
	cerr<<"Program terminated successfully"<<endl;
	return 0;
	break;
    case 'c':  //windows on entire chr
	v=rw.getGenomicWindowsChr(chrToUse,bpToExtract,overlapBetweenWindows);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    case 'o':    //windows on entire genome
	v=rw.getGenomicWindows(bpToExtract,overlapBetweenWindows);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    case 'v':  //entire chr
	v=rw.getChr(chrToUse);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    case 'w':  //entire genome
	v=rw.getGenomeWide();
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    default:
	cerr << "Unknown loci selection:"<< lociSelection <<endl;
	return 1;       
	break;
    }


    for(unsigned int i=0;i<v.size();i++){
	if(!inBedFormat)
	    cout<<v[i]<<endl;
	else
	    cout<<v[i].asBed()<<endl;
    }



    return 0;
}

