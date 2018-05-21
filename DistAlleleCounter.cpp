/*
 * DistAlleleCounter
 * Date: Aug-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "DistAlleleCounter.h"

const string    DistAlleleCounter::dimers[16] = {
    "AA",
    "AC",
    "AG",
    "AT",

    "CA",
    "CC",
    "CG",
    "CT",

    "GA",
    "GC",
    "GG",
    "GT",

    "TA",
    "TC",
    "TG",
    "TT"};

DistAlleleCounter::DistAlleleCounter(){
    reinitializedCounters();



}

void DistAlleleCounter::reinitializedCounters(){
    for(int i=0;i<16;i++)
	count[i] = 0;
}

void DistAlleleCounter::addAllelePair(int indexDimer){
    count[indexDimer]++;
}


string DistAlleleCounter::headerForCount() const{
    string dnaAlphabet = "ACGT";
    string toReturn="";

    for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	    toReturn+=dnaAlphabet.substr(i,1)+dnaAlphabet.substr(j,1);
	    if(j!=3)
		toReturn+="\t";
	}
	if(i!=3)
	    toReturn+="\t";	
    }

    toReturn+="\tident\tmutat\trate";	

    return toReturn;
}

unsigned int DistAlleleCounter::getIdent() const{
    //unsigned int mutations=0;
    unsigned int ident    =0;
    for(int i=0;i<16;i++){
	if( i == 0 || i == 5 || i == 10 || i == 15 ){
	    ident+=     count[i] ;
	}else{
	    //mutations+= apc.count[i];	
	}
    }
    return ident;
}

unsigned int DistAlleleCounter::getMutations() const{
    unsigned int mutations=0;
    for(int i=0;i<16;i++){
	if( i == 0 || i == 5 || i == 10 || i == 15 ){
	    //ident+=     apc.count[i] ;
	}else{
	    mutations+= count[i];	
	}
    }
    return mutations;
}

// double DistAlleleCounter::lowerConf (const unsigned int shortBranch,const unsigned int commonBranch) const{
//     double pHat = double(shortBranch)/ ( double(shortBranch)+double(commonBranch));
//     //1.0 = 85%, 
//     //1.6 = 95%
//     double z = 1.6;
//     double n=shortBranch+commonBranch;
//     //taken from http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
//     double numerator= (  pHat + (z*z/(2.0*n)) - z* sqrt(  (pHat*(1.0-pHat)+(z*z/(4.0*n)) )/n ));
//     double denominator = (1.0+(z*z/n));

//     // cout<<"pHat "<<pHat<<endl;
//     // // cout<<(z*z)<<endl;
//     // // cout<<(2.0*n)<<endl;
//     // cout<<((z*z)/(2.0*n))<<endl;
//     // cout<<numerator<<endl;
//     // cout<<denominator<<endl;

//     // cout<<"####"<<endl;
    
//     return ( numerator /  denominator   );
// }

string DistAlleleCounter::getHeader(string prefixToAdd){
    string dnaAlphabet = "ACGT";
    string toReturn=""+prefixToAdd+"\t";

    for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	    toReturn+=dnaAlphabet.substr(i,1)+dnaAlphabet.substr(j,1);
	    if(j!=3)
		toReturn+="\t";
	}
	if(i!=3)
	    toReturn+="\t";	
    }

    toReturn+="\tident\tmutat\trate";	

    return toReturn;
    // return 
    // 	prefixToAdd+"noMut\t"+
    // 	prefixToAdd+"mut\t";

}

// double DistAlleleCounter::highConf (const unsigned int shortBranch,const unsigned int commonBranch) const {
//     double pHat = double(shortBranch)/ ( double(shortBranch)+double(commonBranch));
//     //1.0 = 85%, 
//     //1.6 = 95%
//     double z = 1.6;
//     double n=shortBranch+commonBranch;
//     //taken from http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
//     double numerator= (  pHat + (z*z/(2.0*n)) + z* sqrt(  (pHat*(1.0-pHat)+(z*z/(4.0*n)) )/n ));
//     double denominator = (1.0+(z*z/n));

//     return ( numerator /  denominator   );
// }

// int main (int argc, char *argv[]) {
//     DistAlleleCounter al;
//     cout<<al.lowerConf(46,396)<<endl;
// }
// double DistAlleleCounter::distRefSam () const {

//     if( (counterReference + counterCommon) != 0){
// 	return double(counterReference)/double(counterReference+counterCommon);
//     }else{
// 	return std::numeric_limits<double>::infinity();
//     }    

// }



// double DistAlleleCounter::distSamRef () const {

//     if( (counterSample + counterCommon) != 0){
// 	return double(counterSample)/double(counterSample+counterCommon);
//     }else{
// 	return std::numeric_limits<double>::infinity();
//     }    

// }
//This subroutine returns as a double the distance as specified by dnaDistMode
double DistAlleleCounter::dist(const string & dnaDistMode) const{
    if(dnaDistMode == "none"){
	double d = double(getMutations())   / double( getIdent()+getMutations() );
	return d;
    }

    if(dnaDistMode == "JC69"){	
	double P = double(getMutations())   / double( getIdent() );
	double inLog=(1.0-P*double(4.0)/double(3.0));
	
	if (inLog != 0.0 )
	    return ( -0.75*log(inLog) );

	return ( 0.0 );
    }


    cerr<<"ERROR: unknown mode for DNA distance: \""<< dnaDistMode<<"\""<<endl;
    exit(1);

    return -1;
}

string DistAlleleCounter::printWithLabel() const{
    stringstream toreturn;

    for(int i=0;i<15;i++){
	toreturn<<dimers[i]<<":"<<count[i]<<"\t";
	//cerr<<dimers[i];
	//toreturn<<count[i]<<"\t";
    }
    toreturn<<dimers[15]<<":"<<count[15];

    return toreturn.str();

}
