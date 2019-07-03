/*
 * SimpleVCF
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "SimpleVCF.h"
//#define DEBUG

//static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

SimpleVCF::SimpleVCF(string line){
    //trimWhiteSpacesBothEnds	(&line);
    vector<string> fields=allTokens(line,'\t');
    corevcf = new CoreVCF(fields);
    deleteCore=true;
    // cerr<<"Ok "<<endl;
    init(fields,corevcf);
    // cerr<<"Ok "<<endl;
    // cout<<"SimpleVCF "<<this<<"\t"<<deleteCore<<endl;
}


SimpleVCF::SimpleVCF(const vector<string> & fields, CoreVCF *  corevcf_,bool deleteCore_):corevcf(corevcf_),deleteCore(deleteCore_){ //string line){
    init(fields,corevcf_);
}

SimpleVCF::SimpleVCF(const char * p, CoreVCF *  corevcf_,bool deleteCore_):corevcf(corevcf_),deleteCore(deleteCore_){ //string line){
    initWithString(p,corevcf_);
}

void SimpleVCF::init_(CoreVCF *  corevcf_){ 
    unresolvedGT=false;
    homozygousREF=false;
    heterozygous=false;
    homozygousALT=false;

    indexGenotype=-1; 
    indexGenotypeQual=-1; 
    indexDepth=-1;    
    indexPL=-1;       

    typeOfData=1;


    // fields=allTokens(line,'\t');
    corevcf = corevcf_;

}

void SimpleVCF::init(const vector<string> & fields, CoreVCF *  corevcf_){ //string line){

    init_(corevcf_);
    
    int fieldIndex  = corevcf->getFieldIndexAndIncrease();
#ifdef DEBUG
    cerr<<"fieldIndex "<<fieldIndex<<endl;
#endif

    //FORMAT FIELDS
    //rawFormatNames  = fields[ corevcf->getFieldIndexINFO()+1 ];
    rawFormatValues = fields[fieldIndex];
    init2();
}

void SimpleVCF::initWithString(const char * p, CoreVCF *  corevcf_){ //string line){

    init_(corevcf_);
    
    //     int fieldIndex  = corevcf->getFieldIndexAndIncrease();
    // #ifdef DEBUG
    //     cerr<<"fieldIndex "<<fieldIndex<<endl;
    // #endif

    //FORMAT FIELDS
    //rawFormatNames  = fields[ corevcf->getFieldIndexINFO()+1 ];
    //rawFormatValues = fields[fieldIndex];
    rawFormatValues = string(p);
    init2();
}

void SimpleVCF::init2(){

#ifdef DEBUG
    //    cerr<<"rawFormatNames  "<<rawFormatNames<<endl;
    cerr<<"rawFormatValues "<<rawFormatValues<<endl;
#endif

    //cout<<rawFormatValues<<endl;
    //formatFieldNames  = allTokens(rawFormatNames ,':');
    //formatFieldNames  = allTokens(rawFormatNames ,':');
    formatFieldNames  = corevcf->getFormatNames();
    formatFieldValues = allTokens(rawFormatValues,':');
    
    if(rawFormatValues == "./."){
	unresolvedGT=true; 
	observedPL=false;
	observedGL=false;
	haploidCall=false;
    }else{

	if(formatFieldNames->size() != formatFieldValues.size()){
	    //cerr<<"SimpleVCF: for line "<<vectorToString(formatFieldValues,":")<<" the format field does not have as many fields as the values "<<formatFieldNames->size()<<" vs "<<formatFieldValues.size() <<endl;
	    //exit(1);
	}

    observedPL=false;
    observedGL=false;
    haploidCall=false;
    for(unsigned int i=0;i<formatFieldNames->size();i++){
      //cerr<<"formatFieldNames["<<i<<"] "<<formatFieldNames->at(i)<<" = "<<formatFieldValues[i]<<endl;
	
	if(formatFieldNames->at(i) == "GT"){ 
	    indexGenotype     =i; 
	    formatFieldGT=                   formatFieldValues[i]; 
	    bool determinedGenotype=false;
	    //Taken from http://www.broadinstitute.org/gatk/guide/topic?name=intro


	    if(formatFieldGT == "./."){ determinedGenotype=true; unresolvedGT=true;	    } 

	    if(formatFieldGT == "0"){   determinedGenotype=true; homozygousREF=true;   haploidCall=true;   }
	    if(formatFieldGT == "1"){   determinedGenotype=true; homozygousALT=true;   haploidCall=true;   }

	    if(formatFieldGT == "0/0"){ determinedGenotype=true; homozygousREF=true;      }
	    if(formatFieldGT == "0|0"){ determinedGenotype=true; homozygousREF=true;      }

	    if(formatFieldGT == "0/1"){ determinedGenotype=true; heterozygous=true;       }
	    if(formatFieldGT == "0|1"){ determinedGenotype=true; heterozygous=true;       }
	    if(formatFieldGT == "1|0"){ determinedGenotype=true; heterozygous=true;       }

	    if(formatFieldGT == "1/1"){ determinedGenotype=true; homozygousALT=true;      }
	    if(formatFieldGT == "1|1"){ determinedGenotype=true; homozygousALT=true;      }

	    if(formatFieldGT == "1/2"){ determinedGenotype=true; heterozygousALT=true;    } //has first alt and second alt
	    if(formatFieldGT == "1|2"){ determinedGenotype=true; heterozygousALT=true;    } //has first alt and second alt

	    if(formatFieldGT == "2/1"){ determinedGenotype=true; heterozygousALT=true;    } //has first alt and second alt
	    if(formatFieldGT == "2|1"){ determinedGenotype=true; heterozygousALT=true;    } //has first alt and second alt

	    if(formatFieldGT == "0|2"){ determinedGenotype=true; heterozygous2ndALT=true; } //has ref       and second alt
	    if(formatFieldGT == "0/2"){ determinedGenotype=true; heterozygous2ndALT=true; } //has ref       and second alt

	    if(formatFieldGT == "2/0"){ determinedGenotype=true; heterozygous2ndALT=true; } //has ref       and second alt
	    if(formatFieldGT == "2|0"){ determinedGenotype=true; heterozygous2ndALT=true; } //has ref       and second alt

	    if(formatFieldGT == "2/2"){ determinedGenotype=true; homozygous2ndALT=true; }   //twice the second alt
	    if(formatFieldGT == "2|2"){ determinedGenotype=true; homozygous2ndALT=true; }   //twice the second alt

	    //for more than 3

	    if(!determinedGenotype){

		vector<string> fieldsOfGT  = allTokens(formatFieldGT ,'/');
	    
		if(fieldsOfGT.size() == 2){
		    // int alleleCFirst = destringify<int>  (fieldsOfGT[0]);
		    // int alleleC2nd   = destringify<int>  (fieldsOfGT[1]);
		    if(isPositiveInt(fieldsOfGT[0]) &&
		       isPositiveInt(fieldsOfGT[1])    ){
			determinedGenotype=true; unresolvedGT=true; 
		    }		   
		}else{
		    vector<string> fieldsOfGT  = allTokens(formatFieldGT ,'|');
	    
		    if(fieldsOfGT.size() == 2){
			// int alleleCFirst = destringify<int>  (fieldsOfGT[0]);
			// int alleleC2nd   = destringify<int>  (fieldsOfGT[1]);
			if(isPositiveInt(fieldsOfGT[0]) &&
			   isPositiveInt(fieldsOfGT[1])   ){
			    determinedGenotype=true; unresolvedGT=true; 
			}		   
		    }else{

		    }

		}
	    }
	    // if(formatFieldGT == "0/3" ||
	    //    formatFieldGT == "3/3" ||
	    //    formatFieldGT == "3/3" ||
	       
	       
	    //    ){ determinedGenotype=true; unresolvedGT=true; 
	   

	    
	    if(!determinedGenotype){
		//cerr<<"SimpleVCF: unable to determine genotype for line "<<vectorToString(fields,"\t")<<" field=#"<<formatFieldGT<<"#"<<endl;
	      cerr<<"SimpleVCF: unable to determine genotype for field=#"<<formatFieldGT<<"#"<<" at position: "<<corevcf->getChr()<<":"<<corevcf->getPosition()<<endl;
	      exit(1);
	    }

	    if(formatFieldNames->size() != formatFieldValues.size()){
		cerr<<"SimpleVCF: WARNING: for the record "<<vectorToString(formatFieldValues,":")<<", the format field does not have as many fields as the values "<<formatFieldNames->size()<<" vs "<<formatFieldValues.size() <<", we used "<<formatFieldGT<<" as genotype and ignored the rest"<<endl;	
		break;
	    }
	    
	    if(formatFieldValues.size() == 1 ){//weird case where the remaining fields are missing
	      break;
	    }
	    continue;
	}



	if(formatFieldNames->at(i) == "GQ"){ 
	    if(formatFieldValues[i] == "."){
		indexGenotypeQual =i; 
		formatFieldGQ=0.0;
	    }else{
		indexGenotypeQual =i; 
		formatFieldGQ=destringify<float>(formatFieldValues[i]);
	    } 
	    continue; }

	if(formatFieldNames->at(i) == "DP"){ 
	    indexDepth        =i;

	    if(formatFieldValues[i] == "."){
		formatFieldDP=-1;//if the DP field is missing
	    }else{
		if(!formatFieldValues[i].empty())
		    formatFieldDP=destringify<int>  (formatFieldValues[i]); 
		else
		    formatFieldDP=-1;//if the DP field is missing
	    }
	    continue;
	}

	if(formatFieldNames->at(i) == "GL"){ 
	    observedGL=true;
	    if(observedPL){
		//cerr<<"SimpleVCF: cannot observed both GL and PL "<<vectorToString(fields,"\t")<<""<<endl;
		cerr<<"SimpleVCF: cannot observed both GL and PL"<<endl;
		exit(1);
	    }

	    indexPL        = i; 
	    formatFieldGL  = formatFieldValues[i];
	    vector<string> glfields = allTokens(formatFieldGL,',');

	    if(glfields.size() == 2){ //haploid calls (e.g. X for a male)
		if(!haploidCall){
		    //cerr<<"SimpleVCF: cannot observed 2 GL fields for a non-haploid record "<<vectorToString(fields,"\t")<<""<<endl;
		    cerr<<"SimpleVCF: cannot observed 2 GL fields for a non-haploid record "<<endl;
		    exit(1);
		}
		formatFieldPLHomoRef =  int(-10.0*destringify<double>(glfields[0]));
		formatFieldPLHetero  =  -1000000; //very unlikely
		formatFieldPLHomoAlt =  int(-10.0*destringify<double>(glfields[1]));
		    
	    }else{
		if(glfields.size() == 3){ //biallelic

		    formatFieldPLHomoRef =  int(-10.0*destringify<double>(glfields[0]));
		    formatFieldPLHetero  =  int(-10.0*destringify<double>(glfields[1]));
		    formatFieldPLHomoAlt =  int(-10.0*destringify<double>(glfields[2]));

		}else{
		    if(glfields.size() == 6){ //triallelic
			//according to VCF docs it has the following order AA,AB,BB,AC,BC,CC
			formatFieldPLHomoRef  =  int(-10.0*destringify<double>(glfields[0])); //r-r

			formatFieldPLHetero1  =  int(-10.0*destringify<double>(glfields[1])); //r-a1
			formatFieldPLHomoAlt1 =  int(-10.0*destringify<double>(glfields[2])); //a1-a1

			formatFieldPLHetero2  =  int(-10.0*destringify<double>(glfields[3])); //r-a2

			formatFieldPLHetero12 =  int(-10.0*destringify<double>(glfields[4])); //a1-a2
			formatFieldPLHomoAlt2 =  int(-10.0*destringify<double>(glfields[5])); //a2-a2

		    }else{
			cerr<<"SimpleVCF: for line "<<vectorToString(glfields,"\t")<<" the GL field does not have 3 or 6 fields"<<endl;
			exit(1);
		    }
		}
	    }
	}


	if(formatFieldNames->at(i) == "PL"){ 
	    observedPL=true;

	    if(observedGL){
		//cerr<<"SimpleVCF: cannot observed both GL and PL "<<vectorToString(fields,"\t")<<""<<endl;
		cerr<<"SimpleVCF: cannot observed both GL and PL"<<endl;
		exit(1);
	    }

	    indexPL        = i; 

	    if(formatFieldValues[i] == "."){
		formatFieldPL = formatFieldValues[i];
		unresolvedGT=true; 
		continue;
	    }

	    formatFieldPL  = formatFieldValues[i];

	    // bizarre format 
	    // the PL are missing if we have a homo ref site.
	    // probably will cause some reference bias
	    // setting the PL to homo ref with infinite quality
	    // when there is no alt with samtools, no PL field are reported
	    if( (formatFieldPL == "0" &&
		 formatFieldGT == "0/0" )
		||
		(corevcf->getAlt()  == ".")
	    ){ 
#ifdef DEBUG
		//		cerr<<"rawFormatNames  "<<rawFormatNames<<endl;
		cerr<<"rawFormatValues "<<rawFormatValues<<endl;
		cerr<<"unresolvedGT    "<<unresolvedGT<<endl;
#endif
		// exit(1);
		formatFieldPLHomoRef =  0;
		formatFieldPLHetero  =  32767;
		formatFieldPLHomoAlt =  32767;

		// pair<int,int> tempp = this->returnLikelyAlleleCountForRefAlt(30);
		// cerr<<tempp.first<<endl;
		// cerr<<tempp.second<<endl;

		continue;
	    }
	       
	    vector<string> plfields = allTokens(formatFieldPL,',');

	    if(plfields.size() == 3){ //biallelic
		formatFieldPLHomoRef =  plString2Int(plfields[0]);
		formatFieldPLHetero  =  plString2Int(plfields[1]);
		formatFieldPLHomoAlt =  plString2Int(plfields[2]);

	    }else{
		if(plfields.size() == 6){ //triallelic
		    //according to VCF docs it has the following order AA,AB,BB,AC,BC,CC
		    formatFieldPLHomoRef  =  plString2Int(plfields[0]); //r-r

		    formatFieldPLHetero1  =  plString2Int(plfields[1]); //r-a1
		    formatFieldPLHomoAlt1 =  plString2Int(plfields[2]); //a1-a1

		    formatFieldPLHetero2  =  plString2Int(plfields[3]); //r-a2

		    formatFieldPLHetero12 =  plString2Int(plfields[4]); //a1-a2
		    formatFieldPLHomoAlt2 =  plString2Int(plfields[5]); //a2-a2

		}else{
		    if(plfields.size() == 10){ //penta allelic
			//according to VCF docs it has the following order  AA,AB,BB,AC,BC,CC, AD,BD,CD,DD
			formatFieldPLHomoRef  =  plString2Int(plfields[0]); //r-r
			
			formatFieldPLHetero1  =  plString2Int(plfields[1]); //r-a1
			formatFieldPLHomoAlt1 =  plString2Int(plfields[2]); //a1-a1

			formatFieldPLHetero2  =  plString2Int(plfields[3]); //r-a2

			formatFieldPLHetero12 =  plString2Int(plfields[4]); //a1-a2
			formatFieldPLHomoAlt2 =  plString2Int(plfields[5]); //a2-a2

			formatFieldPLHetero3  =  plString2Int(plfields[6]); //r-a3
			
			formatFieldPLHetero13 =  plString2Int(plfields[7]); //a1-a3
			formatFieldPLHetero23 =  plString2Int(plfields[8]); //a2-a3
			
			formatFieldPLHomoAlt3 =  plString2Int(plfields[9]); //a1-a1

		    }else{
			//to take care of such things:
			//chr17	55936366	.	N	A,C,T,G	19	q:30	DP=12;VDB=8.329641e-02;AF1=1;AC1=2;DP4=0,0,1,6;MQ=36;FQ=-36	GT:PL:GQ	1/1:63,21,12,51,0,45,58,10,42,55,58,10,42,50,55:15
			if(plfields.size() == 15){
			    unresolvedGT=true; 
			    observedPL=false;
			    observedGL=false;
			    haploidCall=false;
			    return;
			}

			cerr<<"SimpleVCF: for line = "<< corevcf->getChr()<<" "<<corevcf->getPosition()<<" "<<vectorToString(plfields,"\t")<<" the PL field does not have 3,6 or 10 fields, has "<<plfields.size()<<" fields"<<endl;
			exit(1);
		    }
		}
	    }
	    continue;
	}

	//To uncomment the fields to get these fields
	if(formatFieldNames->at(i) == "A"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(unsigned int j=0;j<adfield.size();j++){
		countA.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}

	if(formatFieldNames->at(i) == "C"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(unsigned int j=0;j<adfield.size();j++){
		countC.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}

	if(formatFieldNames->at(i) == "G"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(unsigned int j=0;j<adfield.size();j++){
		countG.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}

	if(formatFieldNames->at(i) == "T"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(unsigned int j=0;j<adfield.size();j++){
		countT.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}
	    


    }
    }
    // cout<<getADforA()<<endl;
    // cout<<getADforC()<<endl;
    // cout<<getADforG()<<endl;
    // cout<<getADforT()<<endl;

    // cerr<<"end"<<endl;

}

SimpleVCF::~SimpleVCF(){

    ///cout<<"DESTRUCTOR SimpleVCF"<<endl;
    if(deleteCore){
	// cout<<"delete CORE "<<this<<endl;
	delete corevcf;
    }
    // delete  altAlleles;
    // delete  fields;
    // delete  infoField;
    // delete formatFieldNames;
    // delete formatFieldValues;
    //exit(1);
    // if(haveInfoField)
    // 	delete infoField;
}



CoreVCF * SimpleVCF::getCorevcf(){
    return corevcf;
}



string  SimpleVCF::getGenotype() const{
    if(indexGenotype != -1){
	return formatFieldGT;
    }else{
	return "";
    }
}


float   SimpleVCF::getGenotypeQual() const{
    if(indexGenotypeQual != -1){
	return formatFieldGQ;
    }else{
	return -1.0;
    }
}

int     SimpleVCF::getDepth() const{
    if(indexDepth != -1){
	return formatFieldDP;
    }else{
	return -1;
    }
}

string  SimpleVCF::getPL() const{
    if(indexPL != -1){
	return formatFieldPL;
    }else{
	return "";
    }
}

int     SimpleVCF::getPLHomoRef() const{
    if(indexPL != -1){
	return formatFieldPLHomoRef;
    }else{
	return -1;
    }
}

int     SimpleVCF::getPLHetero() const{
    if(indexPL != -1){
	return formatFieldPLHetero;
    }else{
	return -1;
    }
}

int     SimpleVCF::getPLHomoAlt() const{
    if(indexPL != -1){
	return formatFieldPLHomoAlt;
    }else{
	return -1;
    }
}


bool SimpleVCF::getObservedGL() const{
    return observedGL;
}

bool SimpleVCF::getObservedPL() const{
    return observedPL;
}


bool SimpleVCF::containsIndel() const{
    return corevcf->isIndel;
}


bool    SimpleVCF::isResolvedSingleBasePairREF() const{
    return corevcf->resolvedSingleBasePairREF;
}

bool    SimpleVCF::isResolvedSingleBasePairALT() const{
    return corevcf->resolvedSingleBasePairALT;
}


bool SimpleVCF::isUnresolvedGT() const{
    return unresolvedGT;
}

bool SimpleVCF::isHomozygousREF() const{
    return homozygousREF;
}

bool SimpleVCF::isHeterozygous() const{
    return heterozygous;
}

bool SimpleVCF::isHomozygousALT() const{
    return homozygousALT;
}

bool SimpleVCF::isHeterozygousALT() const{
    return heterozygousALT;
}

bool SimpleVCF::hasInfoField(string tag) { return corevcf->hasInfoField(tag); } 
int SimpleVCF::getDepthInfo()  { return corevcf->getDepthInfo(); } 
string SimpleVCF::getRef() const { return corevcf->getRef(); } 
string SimpleVCF::getAlt() const { return corevcf->getAlt(); } 
string SimpleVCF::getID() const { return corevcf->getID(); } 
string SimpleVCF::getChr() const { return corevcf->getChr(); } 
unsigned int SimpleVCF::getPosition() const { return corevcf->getPosition(); } 
string SimpleVCF::getFilter() const { return corevcf->getFilter(); } 
string SimpleVCF::getInfoFieldRaw() const { return corevcf->getInfoFieldRaw(); } 
float SimpleVCF::getQual() const { return corevcf->getQual(); } 
void SimpleVCF::setCloseIndel(bool closeIndel) { corevcf->setCloseIndel(closeIndel); } 
bool SimpleVCF::getCloseIndel() const { return corevcf->getCloseIndel(); } 

bool SimpleVCF::isHeterozygous2ndALT() const{
    return heterozygous2ndALT;
}

unsigned int SimpleVCF::getAltAlleleCount() const{
    return corevcf->altAlleles.size();
}

char SimpleVCF::getOtherAllele(char baseToAvoid) const{
    if( heterozygous ){ 
	
	if(corevcf->ref[0]     == baseToAvoid){
	    return     corevcf->alt[0]; 
	}else{
	    if(corevcf->alt[0] == baseToAvoid){
		return corevcf->ref[0]; 
	    }else{
		cerr<<"SimpleVCF: getOtherAllele Cannot get another allele besides "<<baseToAvoid<<" for "<<(*this)<<endl;
		exit(1);    
	    }
	}

    }

    cerr<<"SimpleVCF: getOtherAllele Cannot get another allele besides "<<baseToAvoid<<" for "<<(*this)<<endl;
    exit(1);    
    return 'N';
}

char SimpleVCF::getRandomAllele() const{
    if( homozygousREF ){ return corevcf->ref[0]; }
    if( homozygousALT ){ return corevcf->alt[0]; }
    if( heterozygous ){ 
	//pick an allele at random
	if(randomBool())
	    return corevcf->ref[0]; 
	else
	    return corevcf->alt[0]; 
    }
    //error
    //should we allow heterozygous alt ?
    cerr<<"SimpleVCF: Cannot generate a random allele for "<<(*this)<<endl;
    exit(1);    
}


// ostream& operator<<(ostream& os, const SimpleVCF& smvcf){
void SimpleVCF::print(ostream& os) const{
    // cerr<<"Error"<<endl;
    // exit(1);
    // os<<smvcf.chrName<<"\t"
    //   <<smvcf.position<<"\t"
    //   <<smvcf.ref<<"\t"
    //   <<smvcf.alt<<"\t"
    //   <<smvcf.qual<<"\t"
    //   <<smvcf.filter<<"\t"
    //   <<smvcf.infoFieldRaw<<"\t"
    //   <<smvcf.rawFormatNames<<"\t"
    //   <<smvcf.rawFormatValues;

    os<<corevcf->chrName<<"\t"
      <<corevcf->position<<"\t"
      <<corevcf->ref<<"\t"
      <<corevcf->alt<<"\t"
      <<corevcf->qual<<"\t"
      <<corevcf->filter<<"\t"
      <<corevcf->infoFieldRaw<<"\t"
	//      <<rawFormatNames<<"\t"
      <<rawFormatValues;


    // return os;
}


bool SimpleVCF::isThisAllelePresent(char bp) const {

    //cout<<"r "<<corevcf->resolvedSingleBasePairREF<<" a "<<corevcf->resolvedSingleBasePairALT<<" "<<bp<<" "<<homozygousREF<<" "<<heterozygous<<" "<<homozygousALT<<endl;

    if(corevcf->resolvedSingleBasePairREF && corevcf->resolvedSingleBasePairALT){ //only look at sites with a single bp
	if( homozygousREF ){ return (corevcf->ref[0] == bp); }
	if( homozygousALT ){ return (corevcf->alt[0] == bp); }
	if( heterozygous ){  return ( (corevcf->ref[0] == bp) || (corevcf->alt[0] == bp)); }
    } 
    return false;

}


bool SimpleVCF::hasAtLeastOneA() const  {   
    return isThisAllelePresent('A');
}

bool SimpleVCF::hasAtLeastOneC() const  {
    return isThisAllelePresent('C');
}

bool SimpleVCF::hasAtLeastOneG() const  {
    return isThisAllelePresent('G');
}

bool SimpleVCF::hasAtLeastOneT() const  {
    return isThisAllelePresent('T');
}


bool SimpleVCF::hasAllele(int indexAlle) const  {

    if(indexAlle == 1){return isThisAllelePresent('A');}
    if(indexAlle == 2){return isThisAllelePresent('C');}
    if(indexAlle == 3){return isThisAllelePresent('G');}
    if(indexAlle == 4){return isThisAllelePresent('T');}

    cerr<<"SimpleVCF: hasAllele() request for value: "<<indexAlle<<" cannot be completed, must be between 1 and 4, exiting"<<endl;
    exit(1);

}

string SimpleVCF::getAlleCountBasedOnGT() const{
    string toreturn;
    if(unresolvedGT)
	toreturn="0,0";
    else{
	
	if(homozygousREF)
	    if(haploidCall)		
		toreturn="1,0";
	    else
		toreturn="2,0";
	else{
	    if(heterozygous)
		toreturn="1,1";
	    else{
		if(homozygousALT)
		    if(haploidCall)
			toreturn="0,1";
		    else		
			toreturn="0,2";		
		else{
		    cerr<<"SimpleVCF:getAlleCountBasedOnGT() unresolved genotype"<<endl;
		    exit(1);
		}
	    }
	}
    }
    //if(formatFieldGT == "1|1"){ determinedGenotype=true; homozygousALT=true;   
    toreturn=toreturn+":"+(this->isCpg()?"1":"0");
    return toreturn;
}


pair<int,int> SimpleVCF::returnLikelyAlleleCountForRefAlt(int minPLdiffind) const{
    // if(!observedPL){
    // 	cerr<<"SimpleVCF: returnLikelyAlleleCountForRefAlt() cannot be called is PL value hasn't been defined for "<<*this<<endl;
    // 	exit(1);
    // }

  //cerr<<"SimpleVCF "<<corevcf->position<<" "<<formatFieldPLHomoRef<<" "<<formatFieldPLHetero<<" "<<formatFieldPLHomoAlt<<" "<<minPLdiffind<<" "<<unresolvedGT<<" "<<observedPL<<endl;

    if(unresolvedGT) //unresolved, we cannot infer anything
	return pair<int,int>(0,0);

    if(!observedPL) //unresolved, we cannot infer anything
	return pair<int,int>(0,0);


    if ( (formatFieldPLHetero-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind) {  //high likelihood of homo ref, produce 2 alleles ref
	// refAlleles+=2;
	// altAlleles+=0;
	return pair<int,int>(2,0);
    } else{
	if ((formatFieldPLHetero-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind) {  //high likelihood of homo alt, produce 2 alleles alt
	    // refAlleles+=0;
	    // altAlleles+=2;
	    return pair<int,int>(0,2);
	} else {
	    if ((formatFieldPLHomoRef-formatFieldPLHetero) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHetero) >= minPLdiffind) { //high likelihood of hetero, produce 1 allele of each
		// refAlleles+=1;
		// altAlleles+=1;
		return pair<int,int>(1,1);
	    }else{
		if ((formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoAlt) < minPLdiffind ) { //high likelihood of at least one alt, produce 1 allele alt 
		    // refAlleles+=0;
		    // altAlleles+=1;
		    return pair<int,int>(0,1);
		}else{
		    if ( (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoRef) < minPLdiffind ) { // high likelihood of at least one ref, produce 1 allele ref
			// refAlleles+=1;
			// altAlleles+=0;
			return pair<int,int>(1,0);
		    }else{
			// refAlleles+=0;
			// altAlleles+=0;
			return pair<int,int>(0,0);
		    }
		}
	    }
	}
    }// end all cases
    
}


pair<int,int> SimpleVCF::returnLikelyAlleleCountForRefAltJustGT() const{

    if(unresolvedGT) //unresolved, we cannot infer anything
	return pair<int,int>(0,0);

    
    if(homozygousREF)
	return pair<int,int>(2,0);
    if(heterozygous)
	return pair<int,int>(1,1);
    if(homozygousALT)
	return pair<int,int>(0,2);


    return pair<int,int>(0,0);//cannot infer anything

}





bool SimpleVCF::isHeterozygousUsingPL(int minPLdiffind) const{

    if(unresolvedGT)
	return false; //unresolved

    if ((formatFieldPLHomoRef-formatFieldPLHetero) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHetero) >= minPLdiffind) { //high likelihood of hetero, produce 1 allele of each
	return true;
    }else{
	return false;
    }// end all cases

}



char SimpleVCF::getRandomAlleleUsingPL(int minPLdiffind) const{

    // if(!observedPL){
    // 	cerr<<"SimpleVCF: getRandomAlleleUsingPL() cannot be called is PL value hasn't been defined for "<<*this<<endl;
    // 	exit(1);
    // }

    if(unresolvedGT)
	return 'X'; //unresolved

    if ( (formatFieldPLHetero-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind) {  //high likelihood of homo ref, produce 2 alleles ref
	// cout<<position<<"\t"<<"2,0"<<endl;
	return corevcf->ref[0];
    } else{
	if ((formatFieldPLHetero-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind) {  //high likelihood of homo alt, produce 2 alleles alt
	    // cout<<position<<"\t"<<"0,2"<<endl;
	    return corevcf->alt[0];
	} else {
	    if ((formatFieldPLHomoRef-formatFieldPLHetero) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHetero) >= minPLdiffind) { //high likelihood of hetero, produce 1 allele of each
		// cout<<position<<"\t"<<"1,1"<<endl;
		if(randomBool())
		    return corevcf->ref[0]; 
		else
		    return corevcf->alt[0]; 

	    }else{
		if ((formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoAlt) <minPLdiffind ) { //high likelihood of at least one alt, produce 1 allele alt 
		    // cout<<position<<"\t"<<"0,1"<<endl;
		    return corevcf->alt[0];
		}else{
		    if ( (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoRef) < minPLdiffind ) { // high likelihood of at least one ref, produce 1 allele ref
			// cout<<position<<"\t"<<"1,0"<<endl;
			return corevcf->ref[0]; 
		    }else{
			return 'X'; //unresolved
		    }
		}
	    }
	}
    }// end all cases
}

int SimpleVCF::getADforA() const{
    int toReturn=0;
    for(unsigned int j=0;j<countA.size();j++){
	toReturn+=countA[j];
    }
    return toReturn;
}
 
int SimpleVCF::getADforC() const{
    int toReturn=0;
    for(unsigned int j=0;j<countC.size();j++){
	toReturn+=countC[j];
    }
    return toReturn;
}

int SimpleVCF::getADforG() const{
    int toReturn=0;
    for(unsigned int j=0;j<countG.size();j++){
	toReturn+=countG[j];
    }
    return toReturn;
}
 

int SimpleVCF::getADforT() const{
    int toReturn=0;
    for(unsigned int j=0;j<countT.size();j++){
	toReturn+=countT[j];
    }
    return toReturn;
}


int SimpleVCF::getADforAllele(int indexAlle) const  {

    if(indexAlle == 1){return getADforA();}
    if(indexAlle == 2){return getADforC();}
    if(indexAlle == 3){return getADforG();}
    if(indexAlle == 4){return getADforT();}

    cerr<<"SimpleVCF: getADforAllele() request for value: "<<indexAlle<<" cannot be completed, must be between 1 and 4, exiting"<<endl;
    exit(1);
}
