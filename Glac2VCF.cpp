
#include "Glac2VCF.h"


using namespace std;


Glac2VCF::Glac2VCF(){

}

Glac2VCF::~Glac2VCF(){

}

string Glac2VCF::usage() const{


    string usage=string("")+"glac2vcf  [options]  <ACF file>\n"+
	"This program takes an ACF/GLF file and prints the segregating sites in Treemix format\n"+
	"\n"+
	"Options:\n"+
	"\t--root\t\tPrint the root  (Default "+boolStringify(printRoot)+" )\n"+
	"\t--anc\t\tPrint the anc   (Default "+boolStringify(printAnc)+" )\n"+
	"\t--one\t\tPrint records with a single base as homozygous   (Default "+boolStringify(singleAlleleAsHomo)+" )\n"+
	"\t     \t\tex: 1,0:0 becomes GT=0/0  \n"+
	"\n";
    	   	
    return usage;
}

int Glac2VCF::run(int argc, char *argv[]){



    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage()<<endl;
	return 1;       
    }

    int lastOpt=1;    
    for(int i=1;i<(argc);i++){ 
	
	if((string(argv[i]) == "-")  ){
	    lastOpt=i;
	    break;          
	}
	
	if(string(argv[i])[0] != '-' ){
	    lastOpt=i;
	    break;
	}       	

        if( string(argv[i]) == "--root"  ){
	    printRoot=true;
            continue;
        }

        if( string(argv[i]) == "--anc"  ){
	    printAnc=true;
            continue;
        }

        if( string(argv[i]) == "--one"  ){
	    singleAlleleAsHomo=true;
            continue;
        }

	cerr<<"Error unknown option "<<argv[i]<<endl;
	exit(1);
    }


    GlacParser gp (argv[lastOpt]);
    //gp.isACFormat()){

    vector<string> toprintPop;

    for(unsigned j=0;j<gp.getPopulationsNames()->size();j++){

	if(j == 0){
	    if(!printRoot)
		continue;
	    toprintPop.push_back(gp.getPopulationsNames()->at(j));
	    continue;
	}
	
	if(j == 1){
	    if(!printAnc)
		continue;
	    toprintPop.push_back(gp.getPopulationsNames()->at(j));
	    continue;
	}
	toprintPop.push_back(gp.getPopulationsNames()->at(j));
    }

    string programLine;
    for(int i=0;i<(argc);i++){ 
        programLine+=(string(argv[i])+" ");
    }


    cout<<"##fileformat=VCFv4.2"<<endl;
    cout<<"##FILTER=<ID=PASS,Description=\"All filters passed\">"<<endl;
    cout<<"##glactoolsVersion="<<returnGitHubVersion(argv[-1],"")<<endl;
    cout<<"##glactoolsCommand="<<programLine<<endl;
    cout<<"##reference=N/A"<<endl;
    for(unsigned int i=0;i<gp.getNumberOfChromosomes();i++){
	cout<<"##contig=<ID="<<gp.getChrKnown()[i]<<",length="<<gp.getChrKnownLength()[i]<<">"<<endl;
    }
    // ##ALT=<ID=*,Description="Represents allele(s) other than observed.">
    // ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
    // ##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
    // ##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
    // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    // ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
    // ##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
    // ##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
    // ##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
    // ##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
    // ##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
    // ##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
    // ##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
    // ##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
    if( gp.isACFormat() ){
	cout<<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"<<endl;
	//cout<<"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">"<<endl;
    }else{	
	cout<<"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">"<<endl;
    }
    // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  inC6.q25.bam
	
    cout<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"<<vectorToString( toprintPop," ")<<endl;
    AlleleRecords * arr;
    unsigned int totalRecords=0;
    // unsigned int keptRecords=0;
    
    
    while(gp.hasData()){
    	//cout<<"data"<<endl;
    	arr = gp.getData();
    	if( totalRecords!=0 &&  (totalRecords%10000000) == 0){
    	    cerr<<"Looked at "<<totalRecords<<" records @ "<<arr->chr<<":"<<arr->coordinate<<endl;
    	}
	cout<<gp.getChromosomeName(arr->chri)<<"\t"<<arr->coordinate<<"\t"<<"."<<"\t"<<arr->ref<<"\t"<<arr->alt<<"\t"<<0<<"\t"<<"."<<"\t"<<".";

    	// if(test->alt == 'N'){
    	//     continue;
    	// }

	


    	// string toprint="";
    	// int counterIndRef=0;
    	// int counterIndAlt=0;
	if( gp.isACFormat() ){
	    cout<<"\tGT";
	    for(unsigned j=0;j<arr->vectorAlleles->size();j++){
		if(j == 0){
		    if(!printRoot)
			continue;
		}
		if(j == 1){
		    if(!printAnc)
			continue;
		}
		cout<<"\t"<<arr->vectorAlleles->at(j).toGT(singleAlleleAsHomo);
	    }
	}else{
	    cout<<"\tPL";
	    for(unsigned j=0;j<arr->vectorGLs->size();j++){
		if(j == 0){
		    if(!printRoot)
			continue;
		}
		if(j == 1){
		    if(!printAnc)
			continue;
		}
		cout<<"\t"<<arr->vectorGLs->at(j).toPL();
	    }
	}

	cout<<endl;
	totalRecords++;
    }

    cerr<<"Program "<<argv[0]<<" wrote "<<totalRecords<<", terminated gracefully"<<endl;

    return 0;
}

