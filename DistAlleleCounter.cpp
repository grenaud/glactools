/*
 * DistAlleleCounter
 * Date: Aug-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "DistAlleleCounter.h"


DistAlleleCounter::DistAlleleCounter(){
    reinitializedCounters();
}

void DistAlleleCounter::reinitializedCounters(){
    counterSame      =0;
    counterCommon    =0;
    counterReference =0;
    counterSample    =0;
}

double DistAlleleCounter::lowerConf (const unsigned int shortBranch,const unsigned int commonBranch) const{
    double pHat = double(shortBranch)/ ( double(shortBranch)+double(commonBranch));
    //1.0 = 85%, 
    //1.6 = 95%
    double z = 1.6;
    double n=shortBranch+commonBranch;
    //taken from http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    double numerator= (  pHat + (z*z/(2.0*n)) - z* sqrt(  (pHat*(1.0-pHat)+(z*z/(4.0*n)) )/n ));
    double denominator = (1.0+(z*z/n));

    // cout<<"pHat "<<pHat<<endl;
    // // cout<<(z*z)<<endl;
    // // cout<<(2.0*n)<<endl;
    // cout<<((z*z)/(2.0*n))<<endl;
    // cout<<numerator<<endl;
    // cout<<denominator<<endl;

    // cout<<"####"<<endl;
    
    return ( numerator /  denominator   );
}

string DistAlleleCounter::getHeader(string prefixToAdd){
    return 
	prefixToAdd+"noMut\t"+
	prefixToAdd+"common\t"+
	prefixToAdd+"ind1Spec\t"+
	prefixToAdd+"ind2Spec\t"+
	prefixToAdd+"divInd1M\t"+
	prefixToAdd+"divInd1L\t"+
	prefixToAdd+"divInd1H\t"+
	prefixToAdd+"divInd2M\t"+
	prefixToAdd+"divInd2L\t"+
	prefixToAdd+"divInd2H";

}

double DistAlleleCounter::highConf (const unsigned int shortBranch,const unsigned int commonBranch) const {
    double pHat = double(shortBranch)/ ( double(shortBranch)+double(commonBranch));
    //1.0 = 85%, 
    //1.6 = 95%
    double z = 1.6;
    double n=shortBranch+commonBranch;
    //taken from http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    double numerator= (  pHat + (z*z/(2.0*n)) + z* sqrt(  (pHat*(1.0-pHat)+(z*z/(4.0*n)) )/n ));
    double denominator = (1.0+(z*z/n));

    return ( numerator /  denominator   );
}

// int main (int argc, char *argv[]) {
//     DistAlleleCounter al;
//     cout<<al.lowerConf(46,396)<<endl;
// }
double DistAlleleCounter::distRefSam () const {

    if( (counterReference + counterCommon) != 0){
	return double(counterReference)/double(counterReference+counterCommon);
    }else{
	return std::numeric_limits<double>::infinity();
    }    

}



double DistAlleleCounter::distSamRef () const {

    if( (counterSample + counterCommon) != 0){
	return double(counterSample)/double(counterSample+counterCommon);
    }else{
	return std::numeric_limits<double>::infinity();
    }    

}
