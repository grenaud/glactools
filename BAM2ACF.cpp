/*
 * BAM2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#include "BAM2ACF.h"



//using namespace BamTools;
using namespace std;

#define MAXCOV 1000
#define MINLENGTHFRAGMENT 0
#define MAXLENGTHFRAGMENT 1000000
#define MAXBASEQUAL             64      // maximal base quality score, greater qual score do not make a lot of sense
#define MAXMAPPINGQUAL          37     // maximal mapping quality
#define DEBUGHTS

#define MIN2(a,b) (((a)<(b))?(a):(b))
//#define MIN(a,b) (((a)<(b))?(a):(b))

bool hasLineToPrint;
//string lineToPrint;
AlleleRecords * arToPrint; 
char   offsetQual=33;
static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer
string alphabetHTSLIB = "NACNGNNNTNNNNNNN";

int alphabetHTSLIB2idx [16] = {-1, 0, 1,-1,  // N A C N
			       2,-1,-1,-1,   // G N N N
			       3,-1,-1,-1,   // T N N N
			       -1,-1,-1,-1}; // N N N N

bool onlyPP=false;
int minMQ=0;


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



BAM2ACF::BAM2ACF(){

}

BAM2ACF::~BAM2ACF(){

}



static int read_bamACF(void *data, bam1_t *b){ // read level filters better go here to avoid pileup

    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure

    int ret;
    while (1){
        
	ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
	

	//skip fragments marked as:
	//1: unmapped
	//2: secondary alignment, only retain primary
	//3: QC failed
	//4: PCR duplicates
	//this is not modifiable for now
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

	//if a reads is paired and not properly paired, skip if the user said pp
	if(onlyPP)
	    if( bam_is_paired( b) &&    
		!bam_is_properpair(b) ){ 
		continue; 
	    }
		

	if ( (int)b->core.qual < minMQ ) continue;
        //if ( (int)b->core.qual < aux->min_mapQ ) continue;

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
                         "\t"+"--bed [BED file]"       +"\t" +"Restrict positions to those in a BED file (default: none)\n"+   
                         "\t"+"--pp            "       +"\t" +"Restrict pairs to proper pairs (default: "+booleanAsString(onlyPP)+")\n"+   
                         "\t"+"                "       +"\t" +"ancestral/root alleles for hominin samples\n"+   
	
	
	//"\t"+"--qc"+"\t\t\t"+"Use the reads that have failed quality control (Default: "+boolStringify(useQCFail)+" )\n"+
			"\t"+"--qual"+"\t\t\t"+"Base quality cutoffs on a PHRED scale (Default: "+stringify(minBaseQual)+" )\n"+
			"\t"+"--mq"+"\t\t\t"+"Mapping quality cutoffs on a PHRED scale (Default: "+stringify(minMQ)+" )\n"+
			// "\t\t"+"--het"+"\t"+"Produce two fasta files containing the alleles for het sites (Default: "+boolStringify(printRoot)+" )\n"		

	"\n"+
	"\tBy default, the bam2acf filters reads that are:\n"+
	"\t\t1: unmapped\n"+
	"\t\t2: secondary alignment, only retain primary\n"+
	"\t\t3: QC failed\n"+
	"\t\t4: PCR duplicates\n"+
	"\n";
    return usage;
}
  



int BAM2ACF::run(int argc, char *argv[]){
    hasLineToPrint=false;
    //lineToPrint ="";
    arToPrint= NULL;
    epoFileB=false;

    int lastOpt=1;
    string bedfilename;    
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

	// if( string(argv[i]) == "--qc"){
	//     useQCFail=true;
	//     continue;
	// }

	if( string(argv[i]) == "--pp"){
	    onlyPP=true;
	    continue;
	}


	if( string(argv[i]) == "--bed"){
	    bedfilename = string(argv[i+1]);
	    continue;
	}


	if( string(argv[i]) == "--mq"){
	    minMQ=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if( string(argv[i]) == "--qual"){
	    minBaseQual=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}


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



    //    BamReader reader;

    int i; 
    int n;
    int tid;
    //    int reg_tid;
    int beg;
    int end;
    int pos;
    int *n_plp;
    int baseQ = 0;
    int mapQ = 0;
    int min_len = MINLENGTHFRAGMENT;

    int all = 1; //status = EXIT_SUCCESS, 

    int max_depth = MAXCOV;
    const bam_pileup1_t **plp;
    bool reg = false;//!region.empty();

    void *bed = 0; // BED data structure
    if(!bedfilename.empty()){
    	bed = bed_read(bedfilename.c_str()); // BED or position list file can be parsed now
    }

    bam_hdr_t * headerBAM = NULL; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int last_pos = -1, last_tid = -1, ret;


    data = (aux_t **)calloc(1, sizeof(aux_t*)); // data[i] for the i-th input
    //reg_tid = 0; 
    beg = 0; 
    end = INT_MAX;  // set the default region

    int rf;
    data[0] = (aux_t *)calloc(1, sizeof(aux_t));

    data[0]->fp = sam_open_format(bamFileToOpen.c_str(), "r", NULL); // open BAM

    if (data[0]->fp == NULL) {
	cerr<<"ERROR: Could not open BAM file "<<bamFileToOpen<<""<<endl;
	exit(1);

    }
    rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
    if (baseQ) rf |= SAM_QUAL;

    data[0]->min_mapQ = mapQ;                    // set the mapQ filter
    data[0]->min_len  = min_len;                 // set the qlen filter
    data[0]->hdr = sam_hdr_read(data[0]->fp);    // read the BAM header
        
    if (data[0]->hdr == NULL) {
	cerr<<"ERROR: Could not read header for bamfile "<<bamFileToOpen<<""<<endl;
	exit(1);
    }

    headerBAM = data[0]->hdr;


    string fastaIndex=fastaFile+".fai";

    if( !isFile(fastaFile) ){
	cerr<<"The fasta file "<<fastaFile<<" does not exists"<<endl;
	return 1;	
    }

    if( !isFile(fastaIndex) ){
	cerr<<"The fasta file "<<fastaFile<<"  does not have an index: "<<fastaIndex<<endl;
	return 1;	
    }




    BamTools::Fasta * fastaReference=new BamTools::Fasta();
    if ( !fastaReference->Open(fastaFile , fastaIndex) ){	 
	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<""<<endl;
	 return 1;
     }


     vector<chrinfo> chrFound;
     uint64_t genomeLength;
     readFastaIndex(fastaIndex,chrFound,genomeLength);


    // vector<chrinfo> chrFound;
    //uint64_t genomeLength;
    //  for(int ns=0;ns<faidx_nseq(fai);ns++){
    // 	chrinfo toadd;
    // 	toadd.name    = string( faidx_iseq(fai,ns ) );
    // 	toadd.length  = faidx_seq_len(fai, toadd.name.c_str() );

    // }
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



    //EPO stuff
    string epoFile;
    string epoFileidx;
    //bool epoFileB;

    string epoChr;
    unsigned int epoCoord;

    ReadTabix * rtEPO = NULL;
    string lineFromEPO;
    bool lineLeftEPO;
    bool cpgEPO;
    bool firstLine=true;
    char allel_chimp;
    char allel_anc;

    unsigned int previousPosAlign=0;
    bool hasCpreviousPos=false;



#ifdef AROUNDINDELS
    int prevPos=-1000;
#endif

    
    unsigned int totalBasesL=0;//=cv->getTotalBases();
    unsigned int totalSitesL=0; //=cv->getTotalSites();

    n=1;//number of bam files
    


    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bamACF, (void**)data); // initialization
    if (0 < max_depth)
        bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
    else 
	if (!max_depth)
	    bam_mplp_set_maxcnt(mplp,INT_MAX);
    n_plp = (int *)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp   = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)






    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position	

	unsigned int posAlign    = pos+1;
	int          posAlignInt = int(pos+1);

	//cerr<<endl<<(posAlign)<<" "<<"TNAME "<< headerBAM->target_name[tid]<<" -----------------------"<<endl;
	
	if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= headerBAM->n_targets) continue;     // diff number of @SQ lines per file?
	
	char refC='N';
	char altBase       = 'N';
	
	if ( !fastaReference->GetBase(tid, pos, refC ) ) {
	    cerr << "ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
	    exit(1);
	}

	if(refC == 'N'){//skip unresolved
	    continue;
	}


#ifdef DEBUGHTS
	cerr<<endl<<(posAlign)<<" "<<refC<<" -----------------------"<<endl;
#endif

	//EPO


	if(firstLine){
	    firstLine=false;
	    if(epoFileB){
		rtEPO = new ReadTabix( epoFile.c_str()  , 
				       epoFileidx.c_str()  , 
				       headerBAM->target_name[tid],
				       posAlignInt,
				       INT_MAX ); //the destructor should be called automatically
		BAM2ACF::setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
	    }
	}


	if(epoFileB){
	    
	    if(!lineLeftEPO){
		cerr<<"Error, no data in the EPO file "<< headerBAM->target_name[tid] <<":"<< (posAlign) <<endl;
		exit(1);
	    }
	    
	    if(epoChr != headerBAM->target_name[tid]){
		//cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<m_references[pileupData.RefId].RefName<<endl;
		//exit(1);
		//try to reposition
		rtEPO->repositionIterator( headerBAM->target_name[tid]  , posAlignInt,INT_MAX);
		BAM2ACF::setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
		if(epoChr != headerBAM->target_name[tid]){
		    cerr<<"Error, the repositioning did not work, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<headerBAM->target_name[tid]<<endl;
		    exit(1);
		}
		
	    }
	    


	    while(epoCoord != posAlign){
		if(epoCoord > posAlign){
		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<posAlignInt<<"\t"<<lineFromEPO<<endl;
		    exit(1);
		}
		
		if( ((posAlign) - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO->repositionIterator(headerBAM->target_name[tid]  , pos+1,INT_MAX);
		}

		
		BAM2ACF::setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
	    }
	}
	// #ifdef EPO		    
	// #endif
	//
	// END EPO
	//
	
	int countRef=0;
	int countAlt=0;
	int    basesRetained=0;
	bool foundSites=false;
	//	bool triAllelic=false;
	
        if (all) {
            while (tid > last_tid) {
                if (last_tid >= 0 && !reg) {
                    // Deal with remainder or entirety of last tid.
                    while (++last_pos < int(headerBAM->target_len[last_tid])) {
			//cerr<<(last_pos+1)<<endl;
			//exit(1);
			// // Horribly inefficient, but the bed API is an obfuscated black box.
                        if (bed && bed_overlap(bed, headerBAM->target_name[last_tid], last_pos, last_pos + 1) == 0)
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
                if (bed && bed_overlap(bed, headerBAM->target_name[tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }

            last_tid = tid;
            last_pos = pos;
        }
	if (bed && bed_overlap(bed, headerBAM->target_name[tid], pos, pos + 1) == 0) continue;

	
#ifdef DEBUGHTS
	cerr<<endl<<(posAlign)<<" "<<refC<<" -----------------------"<<endl;
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


        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j;
	    int m = 0;//amount of invalid bases at site j that need to be removed from coverage calculations


#ifdef AROUNDINDELS
	    currIndel   =false;
	    // int  currLevels   [2*MAXCOV];
	    // int  currLevelsi =0;	    
	    currIndelSet.clear();
#endif
	    



            for (j = 0; j < n_plp[i]; ++j) {//each base
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know

#ifdef DEBUGHTS
		cerr<<"i="<<i<<" "<< bam1_qname(p->b)<<" j="<<j<<" "<<int((p->b)->core.flag)<<" "<<p->b->core.n_cigar<<" p="<<(p->cd.p)<<" i="<<int(p->cd.i)<<" f="<<float(p->cd.f)<<" isdel="<<p->is_del<<" isrefsk="<<p->is_refskip<<" indel="<<p->indel<<" "<<p->level<<" "<<p->aux<<" "<<endl;

#endif

		
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


		// This variable will not be zero if the next position is a a deletion of a base in the read or insert in the reference
		// if ( p->indel != 0 ){// having dels or refskips at the next
		//     ++m; 
		//     continue;
		// }

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



		if( j>=MAXCOV){
		    break;
		}

		int bIndex = alphabetHTSLIB2idx[ bam_seqi(bam_get_seq(p->b),p->qpos) ];
		if(bIndex == -1){ continue; } //skip unresolved

		char  b    = "ACGT"[bIndex];
		int   q    = MIN2( int(bam_get_qual(p->b)[p->qpos] ), MAXBASEQUAL);//offsetqual?
		//int   m    = MIN2(  bam_mqual(p->b) , MAXMAPPINGQUAL );
#ifdef DEBUGHTS
		bool isRev = bam_is_rev(p->b);
#endif	
		if(q<minBaseQual){
		    continue;
		}
		//		piToAdd.baseC[bIndex]++;
		//probMM += likeMismatchProbMap[m]; 

		if(b != refC){
		    if(altBase == 'N'){
			altBase = b;
			countAlt++;
		    }else{
			if(b == altBase){
			    countAlt++;
			}else{ //tri-allelic site
			    //triAllelic=true;
			    goto skiptonextpos;
			}
		    }
		}else{
		    countRef++;
		}

		basesRetained++;


		foundSites=true;


#ifdef DEBUGHTS
		//cerr<<isRev<<" "<<"ACGT"[int(sr_.base)]<<" "<<int(sr_.qual)<<" "<<int(sr_.mapq)<<" "<<int(m)<<" "<<int(sr_.lengthF)<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<" "<<p->b->core.n_cigar<<" p="<<(p->cd.p)<<" i="<<int(p->cd.i)<<" f="<<float(p->cd.f)<<" "<<p->indel<<" "<<p->level<<endl;
		//cerr<<"ADD "<<isRev<<" "<<int(sr_.base)<<" "<<int(sr_.qual)<<" "<<int(sr_.mapq)<<" "<<int(m)<<" "<<int(sr_.lengthF)<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<endl;
		cerr<<"ADD "<<isRev<<" "<<bam_seqi(bam_get_seq(p->b),p->qpos)<<" "<<"ACGT"[bIndex]<<" "<<uint8_t(q)<<" "<<int( (p->b)->core.l_qseq )<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<endl;

//<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<endl;
#endif
		
		// if(isRev){
		//     sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].Alignment.Length-pileupData.PileupAlignments[i].PositionInAlignment-1 ); 
		//     sr_.base = 3 - sr_.base;//complement
		// }else{
		//     sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].PositionInAlignment ); 
		// }
		// sr_.isrv=isRev;
		totalBasesL++;
		//piToAdd.readsVec.push_back(sr_);

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
		    cerr<<p->qpos<<"/"<<int(uint8_t( (p->b)->core.l_qseq )) <<",";
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




	    //if(!triAllelic){
	    if(countRef != 0 || countAlt != 0 ){//we add
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
			if(allel_chimp == refC){//no diff between chimp and reference
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
			if(allel_anc == refC){//no diff between ancestor and reference
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
	    
		//TODO put tri-allelic to skip to label
		AlleleRecords * arCurrent = new AlleleRecords (false);
		arCurrent->chri          = chr2index[ headerBAM->target_name[tid] ]; //m_references[pileupData.RefId].RefName];//probably redundant
		arCurrent->coordinate    = posAlign;
		arCurrent->sizePops      = 1;
		arCurrent->ref           = refC;
		arCurrent->alt           = altBase;
		arCurrent->vectorAlleles = new vector<SingleAllele>();
		    
		SingleAllele sample (countRef, countAlt, false);//to set CpG flag
		    
		arCurrent->vectorAlleles->push_back(root);
		arCurrent->vectorAlleles->push_back(anc);
		arCurrent->vectorAlleles->push_back(sample);
		
		if( hasLineToPrint && //case of CpG
		    previousPosAlign == (posAlign-1) &&
		    hasCpreviousPos                  &&
		    ( refC == 'G' || altBase  == 'G')){
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

		
	    }else{ //else	    if(countRef != 0 || countAlt != 0 ){

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
	    hasCpreviousPos  = ( refC == 'C' || altBase  == 'C');


	}//closes pos
	//}//close !triallelic
    
	if( foundSites ){// && !triAllelic ){
	    //piToAdd.avgMQ =  round(-10*log10(probMM/double(basesRetained)));
	    totalSitesL++;
	}

	///piForGenomicWindow->push_back(piToAdd);
	//cerr<<"readsVec size= "<<piToAdd.readsVec.size()<<endl;
#ifdef AROUNDINDELS
	    prevPos = pos;
#endif

    }//end for all pos

        
    //end while mpileup
    
    if (ret < 0){ //status = EXIT_FAILURE;
	cerr<<"Problem parsing region in bamfile "<<bamFileToOpen<<endl;
	exit(1);
    }
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

    // if (all) {
    //     // Handle terminating region
    //     if (last_tid < 0 && reg && all > 1) {
    //         last_tid = reg_tid;
    //         last_pos = beg-1;
    //     }
    //     while (last_tid >= 0 && last_tid < h->n_targets) {
    //         while (++last_pos < int(h->target_len[last_tid])) {
    //             if (last_pos >= end) break;
    //             if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
    //                 continue;
    // 		totalSitesL++;
    //         }
    //         last_tid++;
    //         last_pos = -1;
    //         if (all < 2 || reg)
    //             break;
    //     }
    // }


    //    for (i = 0; i < n && data; ++i) {
    bam_hdr_destroy(data[0]->hdr);
    if (data[0]->fp) sam_close(data[0]->fp);
    hts_itr_destroy(data[0]->iter);
    free(data);
    //}

    
    //free(reg);
    if (bed) bed_destroy(bed);

    
    cerr<<"Program bam2acf terminated gracefully, looked at TODO  reads"<<endl;

     //clean up
     //reader.Close();

     fastaReference->Close();
     delete(gw);
     //     delete cv;
     
     return 0;
}

