/*
 * BAM2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#include "BAM2ACF.h"



using namespace BamTools;
using namespace std;

#define MAXCOV 250
#define MIN(a,b) (((a)<(b))?(a):(b))

bool hasLineToPrint;
//string lineToPrint;
AlleleRecords * arToPrint; 
char   offsetQual=33;
static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer



void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO){

    lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
    if(!lineLeftEPO){
	cerr<<"Error, missing data in the EPO file"<<endl;
	exit(1);
    }

    vector<string> fieldsEPO  = allTokens(lineFromEPO,'\t');
    epoChr                   = fieldsEPO[0];
    epoCoord                 = string2uint(fieldsEPO[1]);					
    if(fieldsEPO[9] == "1")
	cpgEPO=true;		    
    else
	cpgEPO=false;		    



    allel_anc   = fieldsEPO[3][0];//inferred ancestor
    allel_chimp = fieldsEPO[4][0];//chimp;


}


class glacVisitor : public PileupVisitor {
  
    public:
    glacVisitor(const RefVector& references, Fasta * fastaReference,bool useQCFail,const int minBaseQual,const string epoFile_p,bool epoFileB_,GlacWriter * gw_,   const  map<string,uint16_t> * chr2index_)
            : PileupVisitor()
	    , m_references(references)
	    , m_fastaReference(fastaReference)
	    , m_useQCFail(useQCFail)
	    , m_minBaseQual(minBaseQual)
	    , epoFile(epoFile_p)
	    , epoFileB(epoFileB_)
	    , gw(gw_)
	    , chr2index(chr2index_)
	      //, m_out(out)
        { 
	    previousPosAlign=0;
	    hasCpreviousPos = false;

	    epoFileidx= epoFile+".tbi";

	    cpgEPO=false;
	    firstLine=true;
	}
        ~glacVisitor(void) { 
	    
	    if(hasLineToPrint){
		//cout<<lineToPrint<<":0"<<endl;
		//arToPrint
		if(!gw->writeAlleleRecord(arToPrint)){
		    cerr<<"BAM2ACF: error writing header "<<endl;
		    exit(1);
		}
		delete(arToPrint);
	    }
	    hasLineToPrint=true;		

	}
  
    // PileupVisitor interface implementation




    public:
	// prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData) {
	    char referenceBase = 'N';
	    char altBase       = 'N';

    	    unsigned int posAlign    = pileupData.Position+1;
    	    int          posAlignInt = int(pileupData.Position+1);

	    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
		cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
		return;
	    }
	    referenceBase = toupper(referenceBase);

	    //cout<<
	    if(referenceBase == 'N')
		return;


	    //
	    // BEGIN EPO
	    //
	    if(firstLine){
		firstLine=false;
		if(epoFileB){
		    rtEPO = new ReadTabix( epoFile.c_str()  , 
					   epoFileidx.c_str()  , 
					   m_references[pileupData.RefId].RefName,
					   posAlignInt,
					   INT_MAX ); //the destructor should be called automatically
		    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		}
	    }

	    if(epoFileB){
		
		if(!lineLeftEPO){
		    cerr<<"Error, no data in the EPO file "<< m_references[pileupData.RefId].RefName <<":"<< posAlignInt <<endl;
		    exit(1);
		}

		if(epoChr != m_references[pileupData.RefId].RefName){
		    //cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<m_references[pileupData.RefId].RefName<<endl;
		    //exit(1);
		    //try to reposition
		    rtEPO->repositionIterator(m_references[pileupData.RefId].RefName  , posAlignInt,INT_MAX);
		    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		    if(epoChr != m_references[pileupData.RefId].RefName){
			cerr<<"Error, the repositioning did not work, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<m_references[pileupData.RefId].RefName<<endl;
			exit(1);
		    }
		
		}


	    
		while(epoCoord != posAlign){
		    if(epoCoord > posAlign){
			cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(posAlignInt)<<"\t"<<lineFromEPO<<endl;
			exit(1);
		    }

		    if( (posAlign - epoCoord ) >= limitToReOpenFP){ //seeking in the file
			rtEPO->repositionIterator(m_references[pileupData.RefId].RefName  , posAlignInt,INT_MAX);
		    }


		    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		}
	    }
	    //
	    // END EPO
	    //

	    //m_coverageCounter->at(MIN(pileupData.PileupAlignments.size(),MAXCOV))++;
	    int countRef=0;
	    int countAlt=0;
	    string currentLine;
	    string chimpString;
	    string ancString;
		
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
		if(pileupData.PileupAlignments[i].Alignment.IsFailedQC() && 
		   !m_useQCFail){
		    continue;
		}
		//skip deletion in the reads/insertion in the reference
		if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion ){
		    continue;
		}
		char b =       pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		//char q = pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		char q   = int(pileupData.PileupAlignments[i].Alignment.Qualities[pileupData.PileupAlignments[i].PositionInAlignment]-offsetQual);

		if(q<m_minBaseQual)//base on read does not pass quality cutoff
		    continue;

		if(b != referenceBase){
		    if(altBase == 'N'){
			altBase = b;
			countAlt++;
		    }else{
			if(b == altBase){
			    countAlt++;
			}else{ //tri-allelic site
			    goto skiptonextpos;
			}
		    }
		}else{
		    countRef++;
		}

	    }

	    // cerr<<posAlign<<"\talt1\t"<<altBase<<endl;


	    if(countRef != 0 || countAlt != 0 ){
		//		char alt=altBase;//(toprint->getAlt()=="."?'N':toprint->getAlt()[0]);
		SingleAllele root;
		SingleAllele anc;

		
		//unresolved ancestral allele (A,C,G,T)
		if(!epoFileB){ //no epo file
		    // chimpString="0,0:0";
		    // ancString="0,0:0";
		    root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 
		    anc.setRefCount(0);  anc.setAltCount(0);   anc.setIsCpg(false); 

		}else{
		    if(!isResolvedDNA(allel_chimp)){
			//chimpString="0,0:0";
			root.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 

		    }
		    //resolved ancestral allele
		    else{
			if(allel_chimp == referenceBase){//no diff between chimp and reference
			    //chimpString="1,0:"+string(cpgEPO?"1":"0");
			    root.setRefCount(1); root.setAltCount(0);  root.setIsCpg(cpgEPO); 

			}else{
			    if(altBase == 'N'){//no alt defined, the chimp becomes the alt			    
				altBase = allel_chimp;
				//chimpString="0,1:"+string(cpgEPO?"1":"0");
				root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 

			    }else{
				if(altBase == allel_chimp){//alt is chimp
				    //chimpString="0,1:"+string(cpgEPO?"1":"0");
				    root.setRefCount(0); root.setAltCount(1);  root.setIsCpg(cpgEPO); 
				}else{ //tri-allelic site, discard
				    //				continue;
				    goto skiptonextpos;
				}
			    }
			}
		    }
		
		    // cerr<<posAlign<<"\talt2\t"<<altBase<<endl;

		    if(!isResolvedDNA(allel_anc)){
			//ancString="0,0:0";					
			anc.setRefCount(0); root.setAltCount(0);  root.setIsCpg(false); 							    
		    }
		    //resolved ancestral allele
		    else{
			if(allel_anc == referenceBase){//no diff between ancestor and reference
			    //ancString="1,0:"+string(cpgEPO?"1":"0");
			    anc.setRefCount(1); anc.setAltCount(0);  anc.setIsCpg(cpgEPO); 
			}else{
			    if(altBase == 'N'){//no alt defined, the ancestor becomes the alt			    
				altBase = allel_anc;
				//ancString="0,1:"+string(cpgEPO?"1":"0");			    
				anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 

			    }else{
				if(altBase == allel_anc){//alt is ancestor
				    //ancString="0,1:"+string(cpgEPO?"1":"0");
				    anc.setRefCount(0); anc.setAltCount(1);  anc.setIsCpg(cpgEPO); 
				}else{ //tri-allelic site, discard
				    //continue;
				    goto skiptonextpos;
				}
			    }
			}
		    }
		}
		
		// cerr<<posAlign<<"\talt3\t"<<altBase<<endl;

		//currentLine=m_references[pileupData.RefId].RefName + "\t" +stringify(posAlign)+"\t"+ stringify(referenceBase) + ","+stringify(altBase)+"\t"+chimpString+"\t"+ancString+"\t"+stringify(countRef)+","+stringify(countAlt);		
		AlleleRecords * arCurrent = new AlleleRecords (false);
		arCurrent->chri          = chr2index->at(m_references[pileupData.RefId].RefName);//probably redundant
		arCurrent->coordinate    = posAlign;
		arCurrent->sizePops      = 1;
		arCurrent->ref           = referenceBase;
		arCurrent->alt           = altBase;
		arCurrent->vectorAlleles = new vector<SingleAllele>();
		
		SingleAllele sample (countRef, countAlt, false);//to set CpG flag
		
		arCurrent->vectorAlleles->push_back(root);
		arCurrent->vectorAlleles->push_back(anc);
		arCurrent->vectorAlleles->push_back(sample);
		
		// cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
		//     toprint->getRef()<<","<<
		//     alt<<"\t"<<
		//     chimpString<<"\t"<<
		//     ancString<<"\t"<<
		//     pairCount.first<<","<<pairCount.second<<":"<<(toprint->isCpg()?"1":"0")<<endl;	


		if( hasLineToPrint && //case of CpG
		    previousPosAlign == (posAlign-1) &&
		    hasCpreviousPos                  &&
		    ( referenceBase == 'G' || altBase  == 'G')){
		    //currentLine+=":1";
		    //cout<<lineToPrint<<":1"<<endl;
		    //cout<<currentLine<<":1"<<endl;		
		    arToPrint->vectorAlleles->at(2).setIsCpg(true);
		    if(!gw->writeAlleleRecord(arToPrint)){
			cerr<<"BAM2ACF: error writing header "<<endl;
			exit(1);
		    }
		    delete(arToPrint);

		    arCurrent->vectorAlleles->at(2).setIsCpg(true);
		    if(!gw->writeAlleleRecord(arCurrent)){
			cerr<<"BAM2ACF: error writing header "<<endl;
			exit(1);
		    }
		    delete(arCurrent);

		    hasLineToPrint=false;
		    //lineToPrint ="";
		    arToPrint=NULL;
		}else{//no CpG
		    if(hasLineToPrint){
			//cout<<lineToPrint<<":0"<<endl;
			arToPrint->vectorAlleles->at(2).setIsCpg(false);
			if(!gw->writeAlleleRecord(arToPrint)){
			    cerr<<"BAM2ACF: error writing header "<<endl;
			    exit(1);
			}
			delete(arToPrint);

		    }
		    hasLineToPrint=true;		
		    //lineToPrint = currentLine;		
		    arToPrint = arCurrent;
		}

	      
	    }else{

		//no line to print
		if(hasLineToPrint){
		    //cout<<lineToPrint<<":0"<<endl;
		    arToPrint->vectorAlleles->at(2).setIsCpg(false);
		    if(!gw->writeAlleleRecord(arToPrint)){
			cerr<<"BAM2ACF: error writing header "<<endl;
			exit(1);
		    }
		    delete(arToPrint);
		}
		hasLineToPrint=false;		

	    }



	skiptonextpos:
	    previousPosAlign = posAlign;
	    hasCpreviousPos  = ( referenceBase == 'C' || altBase  == 'C');

        }
        
    private:
    RefVector m_references;
    Fasta * m_fastaReference;
    unsigned int previousPosAlign;
    bool hasCpreviousPos;
    bool m_useQCFail;
    int m_minBaseQual;
    

    //EPO stuff
    string epoFile;
    string epoFileidx;
    bool epoFileB;
    GlacWriter * gw;
    string epoChr;
    unsigned int epoCoord;
    const  map<string,uint16_t> * chr2index;
    ReadTabix * rtEPO ;
    string lineFromEPO;
    bool lineLeftEPO;
    bool cpgEPO;
    bool firstLine;
    char allel_chimp;
    char allel_anc;

};

BAM2ACF::BAM2ACF(){

}

BAM2ACF::~BAM2ACF(){

}

string BAM2ACF::usage() const{
    string usage=string("")+" bam2acf [options] <fasta file> <bam file> <name sample> "+
			"\n\nThis program produces a  mistar matrix given a BAM file\n\n"+
	                "\t<fasta file> is the fasta file of the reference genome used for mapping\n"+
	                "\t<name sample> is the name you wish the sample to have\n"+
			"\tOptions\n"+	     
   	                 "\t"+"-u" + "\t\t\t"+"Produce uncopressed output (default: "+booleanAsString(uncompressed)+")\n"+
                         "\t"+"--epo [EPO file]"       +"\t" +"Use file as EPO alignment to set the (default: none)\n"+   
                         "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   
                  

			"\t"+"--qc"+"\t\t\t"+"Use the reads that have failed quality control (Default: "+boolStringify(useQCFail)+" )\n"
			"\t"+"--qual"+"\t\t\t"+"Base quality cutoffs on a PHRED scale (Default: "+stringify(minBaseQual)+" )\n"
			// "\t\t"+"--het"+"\t"+"Produce two fasta files containing the alleles for het sites (Default: "+boolStringify(printRoot)+" )\n"		
	"";
    return usage;
}
  



int BAM2ACF::run(int argc, char *argv[]){
    hasLineToPrint=false;
    //lineToPrint ="";
    arToPrint= NULL;
    int lastOpt=1;

    //all but last arg
    for(int i=1;i<(argc);i++){ 
	
        if((string(argv[i]) == "-")  ){
            lastOpt=i;
            break;          
        }

        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }

        if(strcmp(argv[i],"--epo") == 0 ){
            epoFile=string(argv[i+1]);
            epoFileB=true;
            i++;
            continue;
        }

        if(string(argv[i]) == "-u"){
            uncompressed=true;
            continue;
        }

	if( string(argv[i]) == "--qc"){
	    useQCFail=true;
	    continue;
	}


	if( string(argv[i]) == "--qual"){
	    minBaseQual=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	// if( string(argv[i]) == "--het"){
	//     produceTwo=true;
	// }
	cerr<<"Unknown option "<<string(argv[i])<<endl;
	return 1;
    }



    if(argc < 4 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    if(lastOpt != (argc-3)){
        cerr<<"The last 3 arguments are <fasta file> <bam file> <name sample> "<<endl;
        return 1;               
    }


    string fastaFile    = string(argv[lastOpt+0]);//fasta file
    string bamfiletopen = string(argv[lastOpt+1]);
    string nameSample   = string(argv[lastOpt+2]);//fasta file


    //string epoFile  = string(argv[argc-1]);
    // string epoFileidx = epoFile+".tbi";
    // string epoChr;
    // unsigned int epoCoord;

    // ReadTabix * rtEPO ;
    // string lineFromEPO;
    // bool lineLeftEPO;
    // bool cpgEPO=false;
    // bool firstLine=true;
    // char allel_chimp;
    // char allel_anc;
    // string regionToUse  = string(argv[argc-3]);//chr:st-end


    // size_t dotdotpos = regionToUse.find("-");
    // if(dotdotpos != string::npos ) {
    // 	regionToUse = regionToUse.substr(0,dotdotpos)+".."+regionToUse.substr(dotdotpos+1);
    // }

    BamReader reader;

     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM file"<<bamfiletopen<< endl;
    	return 1;
     }



     // if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
     // 	 cerr << "Could not open input index BAM files." << endl;
     // 	 return 1;
     // }

    // retrieve reference data
     const RefVector  references = reader.GetReferenceData();

     string fastaIndex=fastaFile+".fai";
     Fasta fastaReference;
     if ( !fastaReference.Open(fastaFile , fastaIndex) ){	 
	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<""<<endl;
	 return 1;
     }

     vector<chrinfo> chrFound;
     uint64_t genomeLength;
     readFastaIndex(fastaIndex,chrFound,genomeLength);

    GlacWriter * gw=NULL;

    gw = new GlacWriter(1,    //gp.getSizePops(),
                        false, //gp.isGLFormat(),
                        2,//gp.isACFormat()?2:1,
			1,//compression threads
                        uncompressed);




    string header="";
    //cout<<"#ACF"<<endl;    
    header+="#ACF\n";
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    //cout<<"#PG:"<<programLine<<endl;
    header+="#PG:"+programLine+"\n";;
    header+="#GITVERSION: "+returnGitHubVersion(argv[-1],"")+"\n";;
    //cout<<"#DATE: "<<getDateString()<<endl;
    header+="#DATE: "+getDateString()+"\n";;
    header+="#BAM2ACF:";

    
    map<string,uint16_t> chr2index;
    uint16_t     chrCurrentIndex=0;

    for(unsigned j=0;j<(chrFound.size());j++){
	header+= string("#SQ\tSN:")+chrFound[j].name+"\tLN:"+stringify(chrFound[j].length)+"\n";
        chr2index[chrFound[j].name]=chrCurrentIndex++;
        if(chrCurrentIndex == 0xFFFF){
            cerr<<"Too many chromosomes for this build, more than 65535"<<endl;
            return 1;
        }
    }

    header+="#chr\tcoord\tREF,ALT\troot\tanc\t"+nameSample+"\n";

    if(!gw->writeHeader(header)){
        cerr<<"GlacViewer: error writing header "<<endl;
        exit(1);
    }
    
    glacVisitor* cv = new glacVisitor(references,&fastaReference,useQCFail,minBaseQual,epoFile,epoFileB,gw,&chr2index);
    PileupEngine pileup;
    pileup.AddVisitor(cv);


     BamAlignment al;
     unsigned int numReads=0;
     while ( reader.GetNextAlignment(al) ) {

	 pileup.AddAlignment(al);
	 numReads++;
     }
     
     pileup.Flush();
     
     if( hasLineToPrint ){
	 //cout<<lineToPrint<<":0"<<endl;
	 
	 if(!gw->writeAlleleRecord(arToPrint)){
	     cerr<<"BAM2ACF: error writing header "<<endl;
	     exit(1);
	 }
	 delete(arToPrint);
	 
     }

     cerr<<"Program bam2acf terminated gracefully, looked at "<<numReads<< " reads"<<endl;

     //clean up
     reader.Close();

     fastaReference.Close();
     delete(gw);
     delete cv;
     
     return 0;
}

