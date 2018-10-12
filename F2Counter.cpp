/*
 * F2Counter
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "F2Counter.h"


F2Counter::F2Counter(){
    reinitializedCounters();
}


void F2Counter::reinitializedCounters(){

    f2Sum       =0.0;
    counterSites=0;
    nInds1      =0.0; 
    nInds2      =0.0; 
    hz1         =0.0; 
    hz2         =0.0; 
    nsnp        =0.0; 

}


string F2Counter::headerForCount() const {
    return "f2sum\tcount\tf2";
}

// vector<double> F2Counter::performBoostraps() const{
//     vector<double> boostraps;
//     timeval time;
//     gettimeofday(&time, NULL);

//     srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );
//     //TODO
//     // unsigned int sumADDA = (counterDerAnc+counterAncDer);
//     // for(unsigned int i=0;i<1000;i++){
//     // 	 unsigned int counterAncDerBoot=0;
//     // 	 unsigned int counterDerAncBoot=0;
//     // 	 for(unsigned int j=0;j<(sumADDA);j++){
//     // 	     unsigned int randnum = rand()%sumADDA;
//     // 	     if(randnum<counterAncDer)
//     // 		 counterAncDerBoot++;
//     // 	     else
//     // 		 counterDerAncBoot++;
//     // 	 }
//     // 	 boostraps.push_back( (double(counterAncDerBoot)-double(counterDerAncBoot))/(double(counterAncDerBoot)+double(counterDerAncBoot)) );
//     // }
    
//     return boostraps;
// }


double F2Counter::computeF2() const{
    double f2= (f2Sum)/(double(counterSites));
    return f2;    
}

//update allele counts for the target
void F2Counter::addAlleleCounts(int c1_1,int c2_1,double f_1,int c1_2,int c2_2,double f_2){
    //cerr<<"c1 "<<c1<<" c2 "<<endl;
    //taken from CountData.cpp from treemx

    nInds1   += ((double) c1_1 + (double) c2_1)/2.0;//number of individuals
    nInds2   += ((double) c1_2 + (double) c2_2)/2.0;//number of individuals

    // double f_1 = (double) c1_1 / ( (double) c1_1 + (double) c2_1 );
    // double f_2 = (double) c1_2 / ( (double) c1_2 + (double) c2_2 );

    hz1 +=  (f_1*(1-f_1));
    hz2 +=  (f_2*(1-f_2));
    // cerr<<"f_1 "<<f_1<<" f_2 "<<f_2<<endl;
    // cerr<<"hz1 "<<hz1<<" hz2 "<<hz2<<endl;

    // double tmp2 = (double) c2 / ((double) c1+ (double) c2 - 1.0);
    // double tmphzy = 2* f * tmp2;
    // if((c1+c2) < 2){
    // 	tmphzy = 2*f*(1-f);
    // }
    // ///cerr<<"c1 "<<c1<<" c2 "<<c2<<" "<<tmphzy<<endl;
    // hzy += tmphzy; //2*f*(1-f);
    nsnp++;
}


F2Counter &  F2Counter::operator+=(const F2Counter & other){
    //As f2 is:
    //we should be able to add substract from the f2Sum and counterSites
    this->f2Sum        += other.f2Sum;
    this->counterSites += other.counterSites;
    this->nInds1       += other.nInds1;
    this->nInds2       += other.nInds2;

    this->hz1          += other.hz1;
    this->hz2          += other.hz2;

    this->nsnp         += other.nsnp;


    return *this;
}


F2Counter &  F2Counter::operator-=(const F2Counter & other){


    this->f2Sum        -= other.f2Sum;
    this->counterSites -= other.counterSites;
    this->nInds1       -= other.nInds1;
    this->nInds2       -= other.nInds2;

    this->hz1          -= other.hz1;
    this->hz2          -= other.hz2;

    this->nsnp         -= other.nsnp;

    return *this;
}
