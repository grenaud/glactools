/*
 * BAM2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#include "BAM2ACF.h"



//using namespace BamTools;
using namespace std;

#define MAXCOV 250
#define MINLENGTHFRAGMENT 0
#define MAXLENGTHFRAGMENT 1000000
//#define MIN(a,b) (((a)<(b))?(a):(b))

bool hasLineToPrint;
//string lineToPrint;
AlleleRecords * arToPrint; 
char   offsetQual=33;
static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer




static int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len) {
    mplp_ref_t *r = ma->ref;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!r || !ma->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?

    if (tid == r->ref_id[0]) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }
    if (tid == r->ref_id[1]) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp;
        tmp = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    free(r->ref[1]);
    r->ref[1]     = r->ref[0];
    r->ref_id[1]  = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                                ma->h->target_name[r->ref_id[0]],
                                0,
                                INT_MAX,
                                &r->ref_len[0]);

    if (!r->ref[0]) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

void BAM2ACF::setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO){

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

/*
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
		    BAM2ACF::setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
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
		    BAM2ACF::setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
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


		    BAM2ACF::setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
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
			anc.setRefCount(0); anc.setAltCount(0);  anc.setIsCpg(false); 							    
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
			    cerr<<"BAM2ACF: error writing allele record "<<endl;
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
			cerr<<"BAM2ACF: error writing allele record "<<endl;
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
*/


BAM2ACF::BAM2ACF(){

}

BAM2ACF::~BAM2ACF(){

}

//TODO remove proper pair to add all
static int read_bamACF(void *data, bam1_t *b){ // read level filters better go here to avoid pileup
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1){
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
	// int32_t qlen   = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
	int32_t qlen   = b->core.l_qseq;

	//int32_t isize  = b->core.n_cigar;

	//skip fragments marked as:
	//1: unmapped
	//2: secondary alignment, only retain primary
	//3: QC failed
	//4: PCR duplicates
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

	//if a reads is paired and not properly paired, skip
	if( bam_is_paired( b) &&    
	    !bam_is_properpair(b) ){ 
	    continue; 
	}
		


        if ( (int)b->core.qual < aux->min_mapQ ) continue;

	// cerr<<aux->min_len<<" "<<qlen<<" "<<int(MINLENGTHFRAGMENT)<<" "<<int(MAXLENGTHFRAGMENT)<<endl;
        if ( aux->min_len && ( (qlen <  int(MINLENGTHFRAGMENT) ) || (qlen > int(MAXLENGTHFRAGMENT) ) )) continue;

	// if(specifiedDeam){
	//     if( bam_is_paired( b) ) continue; //skip paired-end reads because we cannot get proper deamination

	//     uint8_t *rgptr = bam_aux_get(b, "RG");
	//     // cerr<<"rg1 "<<rgptr<<endl;
	//     // cout<<"isize "<<isize<<endl;
	//     string rg = "UNKNOWN";
	    
	//     if(rgptr){
	// 	rg = string( (const char*)(rgptr+1));
	//     }
	//     //cerr<<"rg2 "<<rg<<endl;
	//     //TODO
	//     if(rg2info.find(rg) == rg2info.end()){
	// 	cerr << "Heterozygosity computation: found an unknown RG tag from  "<<bam1_qname(b)<<" for BAM file:" << bamFileToOpen <<endl;
	// 	//exit(1);	
	//     }else{
	// 	if(!rg2info[rg].isPe){ //if single end
	// 	    if(qlen == rg2info[rg].maxReadLength){//probably reached the end of the read length, we will keep maxlength-1 and under		
	// 		//cerr<<"name "<<bam1_qname(b)<<"\t"<<rg<<endl;
	// 		continue;
	// 	    }
	// 	}else{
	// 	    //since we skipped the paired reads, if we have reached here
	// 	    //we can add single-end fragments 		    
	// 	}
	//     }	    
	//     // }else{
	//     // 	//cerr<<"WARNING: Could not retrieve RG tag from : "<<bam1_qname(b)<<" discarding"<<endl;	  	    
	//     // }
	// }
	// if(specifiedDeam){
	//     if( bam_is_paired( b) ) continue; //skip paired-end reads because we cannot get proper deamination
	// }


        break;
    }
    return ret;
}

string BAM2ACF::usage() const{
    string usage=string("")+" bam2acf [options] <fasta file> <bam file> <name sample> "+
			"\n\nThis program produces a  ACF matrix given a BAM file\n\n"+
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
    string bamFileToOpen= string(argv[lastOpt+1]);
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

    //    BamReader reader;

    //if ( !reader.Open(bamfiletopen) ) {
    //cerr << "Could not open input BAM file"<<bamfiletopen<< endl;
    // 	return 1;
    //}

    aux_t *data   =NULL;//bam reader
    hts_idx_t *idx=NULL; //bam index

    data = (aux_t *)calloc(1, sizeof(aux_t));


    data->fp = sam_open_format(bamFileToOpen.c_str(), "r", NULL); // open BAM

    if(data->fp == NULL) {
	cerr<<"ERROR: Could not open input BAM file "<<bamFileToOpen<<""<<endl;
	exit(1);
    }

    idx = sam_index_load(data->fp, bamFileToOpen.c_str());  // load the index
    if (idx == NULL) {
	cerr<<"ERROR: Cannot load index for bamfile "<<bamFileToOpen<<""<<endl;
	exit(1);
    }
    data->hdr = sam_hdr_read(data->fp);    // read the BAM header
    if (data->hdr == NULL) {
	cerr<<"ERROR: Could not read header for bamfile "<<bamFileToOpen<<""<<endl;
	exit(1);
    }



     // if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
     // 	 cerr << "Could not open input index BAM files." << endl;
     // 	 return 1;
     // }

    // retrieve reference data
    //const RefVector  references = reader.GetReferenceData();

    string fastaIndex=fastaFile+".fai";

    if( !isFile(fastaFile) ){
	cerr<<"The fasta file "<<fastaFile<<" does not exists"<<endl;
	return 1;	
    }

    if( !isFile(fastaIndex) ){
	cerr<<"The fasta file "<<fastaFile<<"  does not have an index: "<<fastaIndex<<endl;
	return 1;	
    }

     // Fasta fastaReference;
     // if ( !fastaReference.Open(fastaFile , fastaIndex) ){	 
     // 	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<""<<endl;
     // 	 return 1;
     // }


    faidx_t *fai = fai_load(fastaFile.c_str());
    if (fai == NULL) { 
	cerr << "ERROR: failed to open fasta index " << fastaIndex<<""<<endl;
	return 1;
    }

    
    // typedef struct {
    //     int id; // faidx_t->name[id] is for this struct.
    //     uint32_t line_len, line_blen;
    //     uint64_t len;
    //     uint64_t seq_offset;
    //     uint64_t qual_offset;
    // } faidx1_t;


    vector<chrinfo> chrFound;
    //uint64_t genomeLength;
    for(int ns=0;ns<faidx_nseq(fai);ns++){
	chrinfo toadd;
	toadd.name    = string( faidx_iseq(fai,ns ) );
	toadd.length  = faidx_seq_len(fai, toadd.name.c_str() );

    }
     // 
     // readFastaIndex(fastaIndex,chrFound,genomeLength);

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
    header+="#BAM2ACF:\n";

    
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
        return 1;
    }




    int i, n, tid, reg_tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = MINLENGTHFRAGMENT;
    int prevPos=-1000;

    int all = 0; //status = EXIT_SUCCESS, 

    int max_depth = 20*MAXCOV;
    const bam_pileup1_t **plp;
    //bool reg = !region.empty();

    void *bed = 0; // BED data structure
    string bedfilename;
    if(!bedfilename.empty()){
	bed = bed_read(bedfilename.c_str()); // BED or position list file can be parsed now
    }

    bam_hdr_t *h = NULL; // BAM header of the 1st input
    aux_t **dataIth;
    bam_mplp_t mplp;
    int last_pos = -1, last_tid = -1, ret;

    
     unsigned int totalBasesL=0;//=cv->getTotalBases();
     unsigned int totalSitesL=0; //=cv->getTotalSites();

    n =1;
    
    // data = (aux_t **)calloc(1, sizeof(aux_t*)); // data[i] for the i-th input
    aux_t **    dataArray = (aux_t **)calloc(1, sizeof(aux_t*)); // data[i] for the i-th input
    dataArray[0] = data;
    reg_tid = 0; beg = 0; end = INT_MAX;  // set the default region

    // for (i = 0; i < n; ++i) {
    // 	//cerr<<"i ="<<i<<endl;
    // 	//for (i = 0; i < 1; ++i) {
    //     int rf;
    //     data[i] = (aux_t *)calloc(1, sizeof(aux_t));


    // 	data[i]->fp = sam_open_format(bamfilename.c_str(), "r", NULL); // open BAM

    //     if (data[i]->fp == NULL) {
    // 	    cerr<<"ERROR: Could not open BAM file "<<bamfilename<<""<<endl;
    // 	    exit(1);

    //     }
    //     rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
    //     if (baseQ) rf |= SAM_QUAL;

    //     data[i]->min_mapQ = mapQ;                    // set the mapQ filter
    //     data[i]->min_len  = min_len;                 // set the qlen filter
    //     data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        
    // 	if (data[i]->hdr == NULL) {
    // 	    cerr<<"ERROR: Could not read header for bamfile "<<bamfilename<<""<<endl;
    // 	    exit(1);
    //     }
    //     // if (reg) { // if a region is specified
    // 	//     hts_idx_t *idx = sam_index_load(data[i]->fp, bamfilename.c_str());  // load the index
    //     //     if (idx == NULL) {
    //     //         //print_error("depth", "can't load index for \"%s\"", argv[optind+i]);
    // 	// 	cerr<<"ERROR: Cannot load index for bamfile "<<bamfilename<<""<<endl;
    // 	// 	exit(1);
    //     //     }else{
    // 	// 	//cerr<<"index ok"<<endl;
    // 	//     }
    //     //     data[i]->iter = sam_itr_querys(idx, data[i]->hdr, region.c_str()); // set the iterator
    //     //     hts_idx_destroy(idx); // the index is not needed any more; free the memory
    //     //     if (data[i]->iter == NULL) {
    // 	// 	cerr<<"ERROR: Cannot parse region \""<<region<<"\""<<endl;
    // 	// 	exit(1);
    //     //     }else{
    // 	// 	//cerr<<"region ok"<<endl;
    // 	//     }
    //     // }
    // }

    //h = data[0]->hdr; // easy access to the header of the 1st BAM
    h = data->hdr; // easy access to the header of the 1st BAM
    //dataToWrite->refID = bam_get_tid(h,currentChunk->rangeGen.getChrName().c_str());

    /*
    if (reg) {
        // beg     = data[0]->iter->beg; // and to the parsed region coordinates
        // end     = data[0]->iter->end;
        // reg_tid = data[0]->iter->tid;
        beg     = data->iter->beg; // and to the parsed region coordinates
        end     = data->iter->end;
        reg_tid = data->iter->tid;
    }
    */

    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bamACF, (void**)dataArray); // initialization
    if (0 < max_depth)
        bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp,INT_MAX);
    n_plp = (int *)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp   = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    //int seqptr=0;

#ifdef AROUNDINDELS
    bool prevIndel   =false;
    // int  prevLevels   [2*MAXCOV];
    // int  prevLevelsi =0;
    set<string> prevIndelSet;

    bool currIndel   =false;
    // int  currLevels   [2*MAXCOV];
    // int  currLevelsi =0;	    
    set<string> currIndelSet;

	    
#endif
    
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position	
	//cerr<<endl<<(pos+1)<<" "<<" -----------------------"<<endl;
	
	if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= h->n_targets) continue;     // diff number of @SQ lines per file?
	
	char refC; //= toupper(seq[ pos-currentChunk->rangeGen.getStartCoord()+1 ]);
	int ref_len;<
	mplp_get_ref(data, tid, &ref, &ref_len);

	if(refC == 'N'){//skip unresolved
	    continue;
	}

	// #ifdef DEBUGHTS
	// 	cerr<<endl<<(pos+1)<<" "<<refc<<" -----------------------"<<endl;
	// #endif
	
        if (all) {
            while (tid > last_tid) {
                if (last_tid >= 0 && !reg) {
                    // Deal with remainder or entirety of last tid.
                    while (++last_pos < int(h->target_len[last_tid])) {
			//cerr<<(last_pos+1)<<endl;
			//exit(1);
			// // Horribly inefficient, but the bed API is an obfuscated black box.
                        if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                            continue;
			totalSitesL++;
                    }
                }
                last_tid++;
                last_pos = -1;
                if (all < 2)
                    break;
            }

            // Deal with missing portion of current tid
            while (++last_pos < pos) {
                if (last_pos < beg) continue; // out of range; skip
                if (bed && bed_overlap(bed, h->target_name[tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }

            last_tid = tid;
            last_pos = pos;
        }
	if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue;

	
#ifdef DEBUGHTS
	cout<<endl<<(pos+1)<<" "<<refC<<" -----------------------"<<endl;
#endif
	
#ifdef AROUNDINDELS
	//if indels are found, remove 
	if( (prevPos+1) == pos){
	    prevIndel     = currIndel;
	    prevIndelSet  = currIndelSet;
	    // prevLevelsi = currLevelsi;
	    // for(int ci=0;ci<currLevelsi;ci++){
	    //     prevLevels[ci] = currLevels[ci];
	    // }		
	}else{
	    prevIndel   =false;
	    //prevLevelsi =0;
	    prevIndelSet.clear();
	}
#endif

	//totalSitesL++;

	//fputs(h->target_name[tid], stdout); 
	//printf("\t%d\t", pos+1); // a customized printf() would be faster
	//cerr<<(pos+1)<<endl;
	//previous indels

        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j;
	    int m = 0;//amount of invalid bases at site j that need to be removed from coverage calculations

	    
	    positionInformation piToAdd;

	    piToAdd.skipPosition                 = false;
	    piToAdd.posAlign                     = (pos+1);
	    piToAdd.refBase                      = refC;
	    double probMM=0;
	    int    basesRetained=0;
	    bool foundSites=false;

	    for(int n=0;n<4;n++){
		piToAdd.baseC[n]=0;
	    }
	    
#ifdef AROUNDINDELS
	    currIndel   =false;
	    // int  currLevels   [2*MAXCOV];
	    // int  currLevelsi =0;	    
	    currIndelSet.clear();
#endif
	    
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
		// if( (pos+1) == 16050348){
		//     cerr<<bam1_qname(p->b)<<endl;	  
		// }
#ifdef DEBUGHTS
		cerr<< bam1_qname(p->b)<<" "<<j<<" "<<int((p->b)->core.flag)<<" "<<p->b->core.n_cigar<<" p="<<(p->cd.p)<<" i="<<int(p->cd.i)<<" f="<<float(p->cd.f)<<" "<<p->is_del<<" "<<p->is_refskip<<" "<<p->indel<<" "<<p->level<<" "<<p->aux<<" "<<endl;

#endif

#ifdef AROUNDINDELS	    
		//cerr<<"prevIndel1: "<<prevIndel<<" "<<prevIndelSet.size()<<endl;
#endif
		//prevIndel<<endl;
		
		//base level filters go here
		if (p->is_del || p->is_refskip){

#ifdef AROUNDINDELS	    
		    //currLevels[currLevelsi++] = p->level;
		    currIndel   = true;
		    //cerr<<"inserting "<<bam1_qname(p->b)<<endl;
		    currIndelSet.insert( bam1_qname(p->b) );
#endif
		    ++m; // having dels or refskips at tid:pos
		    continue;
		}

		//cerr<<"prevIndel2: "<<prevIndel<<endl;

		if ( p->indel != 0 ){// having dels or refskips at the next
		    ++m; 
		    continue;
		}

#ifdef AROUNDINDELS

		//cerr<<"prevIndel3: "<<prevIndel<<" "<<prevIndelSet.size()<<endl;
		
		if(prevIndel){
		    //cerr<<"prevIndel, testing  "<<bam1_qname(p->b)<<""<<endl;
		    
		    // for(set<string>::iterator it = prevIndelSet.begin(); it != prevIndelSet.end(); it++)  {
		    //  	cerr<<"record: "<<*it<<endl;
		    //  }
		    
		    if( prevIndelSet.find(bam1_qname(p->b)) != prevIndelSet.end() ){//skip reads previously flagged as indels. Since indels are few, this should be a quicker solution than to read ahead 
			//cerr<<"found  "<<bam1_qname(p->b)<<" in a previous indel"<<endl;			
			continue;
		    }
		    // bool skippos=false;
		    // for(int pi=0;pi<prevLevelsi;pi++){
		    // 	if(p->level == prevLevels[prevLevelsi]){
		    // 	    skippos=true;
		    // 	    break;
		    // 	}
		    // }
		    // if(skippos){
		    // 	continue;
		    // }
		}
#endif		



		// else{ 
		//     if(bam_get_qual(p->b)[p->qpos] < baseQ) 
		// 	++m; // low base quality
		// }
		//TODO reenable
		if( j>=MAXCOV){
		    break;
		}

		int bIndex = alphabetHTSLIB2idx[ bam_seqi(bam_get_seq(p->b),p->qpos) ];
		if(bIndex == -1){ continue; } //skip unresolved

		int   q    = MIN2( int(bam_get_qual(p->b)[p->qpos] ), MAXBASEQUAL);//offsetqual?
		int   m    = MIN2(  bam_mqual(p->b) , MAXMAPPINGQUAL );
		bool isRev = bam_is_rev(p->b);

		piToAdd.baseC[bIndex]++;
		probMM += likeMismatchProbMap[m]; 
		basesRetained++;


		foundSites=true;

		singleRead sr_;
		sr_.base    = uint8_t(bIndex);
		sr_.qual    = uint8_t(q);	    
		sr_.mapq    = uint8_t(m);
		sr_.lengthF = uint8_t( (p->b)->core.l_qseq );
		
		if(isRev){
		    sr_.pos5p = uint8_t( sr_.lengthF - p->qpos -1 );
		    sr_.base  = 3 - sr_.base;//complement
		}else{
		    sr_.pos5p = uint8_t( p->qpos ); 
		}
		sr_.isrv=isRev;

#ifdef DEBUGHTS
		//cerr<<isRev<<" "<<"ACGT"[int(sr_.base)]<<" "<<int(sr_.qual)<<" "<<int(sr_.mapq)<<" "<<int(m)<<" "<<int(sr_.lengthF)<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<" "<<p->b->core.n_cigar<<" p="<<(p->cd.p)<<" i="<<int(p->cd.i)<<" f="<<float(p->cd.f)<<" "<<p->indel<<" "<<p->level<<endl;
		cerr<<"ADD "<<isRev<<" "<<int(sr_.base)<<" "<<int(sr_.qual)<<" "<<int(sr_.mapq)<<" "<<int(m)<<" "<<int(sr_.lengthF)<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<endl;
#endif
		
		// if(isRev){
		//     sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].Alignment.Length-pileupData.PileupAlignments[i].PositionInAlignment-1 ); 
		//     sr_.base = 3 - sr_.base;//complement
		// }else{
		//     sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].PositionInAlignment ); 
		// }
		// sr_.isrv=isRev;
		totalBasesL++;
		piToAdd.readsVec.push_back(sr_);

//currLevelVec.push_back(p->level);
		//totalBasesL ++;

		
#ifdef DEBUGHTS
		if(0){
		// printf("%d,",p->b);
		// printf("%d,",bam_get_seq(p->b));
		// @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
		//  8 for T and 15 for N. Two bases are packed in one byte with the base
		//  at the higher 4 bits having smaller coordinate on the read. It is
		//  recommended to use bam_seqi() macro to get the base.
		//  */
		//coord
		//printf("%d,",p->qpos);
    cerr<<p->qpos<<"/"<<int(sr_.lengthF) <<",";
		//flag
		//printf("%d,",(p->b)->core.flag);
		cerr<<((p->b)->core.flag)<<",";
		//base
		//printf("%d,",bam_seqi(bam_get_seq(p->b),p->qpos));
		//printf("%c,",alphabetHTSLIB[ bam_seqi(bam_get_seq(p->b),p->qpos) ]);
		cerr<<alphabetHTSLIB[ bam_seqi(bam_get_seq(p->b),p->qpos) ]<<",";
		
		//int c = seq_nt16_int[bam_seqi(seq, p->qpos)];
		//qual
		//printf("%d,",bam_get_qual(p->b)[p->qpos]);
		cerr<<int(bam_get_qual(p->b)[p->qpos])<<",";
		//strand
		//printf("%d,",bam_is_rev(p->b));
		cerr<<bam_is_rev(p->b)<<",";
		
		//paired
		//printf("%d,",bam_is_paired(p->b));
		cerr<<bam_is_paired(p->b)<<",";
		cerr<<"D="<<p->is_del<<",RS="<<p->is_refskip<<",";
		//fail
		//printf("%d-",bam_is_failed(p->b));
		cerr<<bam_is_failed(p->b)<<",";
		cerr<<bam_mqual(p->b)<<"-";
		}
#endif
            }//end for each base at pos
            //printf("\tc=%d", n_plp[i] - m); // this the depth to output
	    //cout<<(n_plp[i] - m); // this the depth to output
	    // totalSitesL ++ ;

	    
	    if( foundSites ){
		piToAdd.avgMQ =  round(-10*log10(probMM/double(basesRetained)));
		totalSitesL++;
	    }

	    piForGenomicWindow->push_back(piToAdd);
	    //cerr<<"readsVec size= "<<piToAdd.readsVec.size()<<endl;
	    prevPos = pos;

	}//end for all pos


	// #ifdef DEBUGHTS
	// 	cerr<<endl<<" -----------------------"<<endl;
	// #endif
        
        //putchar('\n');
	//cout<<endl;
    }//end while mpileup
    
    if (ret < 0){ //status = EXIT_FAILURE;
	cerr<<"Problem parsing region:"<<region<<" in bamfile "<<bamfilename<<endl;
	exit(1);
    }
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

    if (all) {
        // Handle terminating region
        if (last_tid < 0 && reg && all > 1) {
            last_tid = reg_tid;
            last_pos = beg-1;
        }
        while (last_tid >= 0 && last_tid < h->n_targets) {
            while (++last_pos < int(h->target_len[last_tid])) {
                if (last_pos >= end) break;
                if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }
            last_tid++;
            last_pos = -1;
            if (all < 2 || reg)
                break;
        }
    }

    //depth_end:
    for (i = 0; i < n && data[i]; ++i) {
        bam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); 
    free(seq);
    fai_destroy(fai);

    
    //free(reg);
    if (bed) bed_destroy(bed);

    
    // glacVisitor* cv = new glacVisitor(references,&fastaReference,useQCFail,minBaseQual,epoFile,epoFileB,gw,&chr2index);
    // PileupEngine pileup;
    // pileup.AddVisitor(cv);


    //  BamAlignment al;
    //  unsigned int numReads=0;
    //  while ( reader.GetNextAlignment(al) ) {

    // 	 pileup.AddAlignment(al);
    // 	 numReads++;
    //  }
     
    //  pileup.Flush();
     
    //  if( hasLineToPrint ){
    // 	 //cout<<lineToPrint<<":0"<<endl;
	 
    // 	 if(!gw->writeAlleleRecord(arToPrint)){
    // 	     cerr<<"BAM2ACF: error writing header "<<endl;
    // 	     exit(1);
    // 	 }
    // 	 delete(arToPrint);
	 
    //  }

     cerr<<"Program bam2acf terminated gracefully, looked at "<<numReads<< " reads"<<endl;

     //clean up
     //reader.Close();

     //fastaReference.Close();
     delete(gw);
     delete cv;
     
     return 0;
}

