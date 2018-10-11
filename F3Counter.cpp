/*
 * F3Counter
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "F3Counter.h"


F3Counter::F3Counter(){
    reinitializedCounters();
}


void F3Counter::reinitializedCounters(){

    f3Sum       =0.0;
    counterSites=0;
    nInds       =0.0; 
    hzy         =0.0; 
    nsnp        =0.0; 

}


string F3Counter::headerForCount() const {
    return "f3sum\tcount\tf3";
}

// vector<double> F3Counter::performBoostraps() const{
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


double F3Counter::computeF3() const{
    double f3= (f3Sum)/(double(counterSites));
    return f3;    
}

//update allele counts for the target
void F3Counter::addAlleleCounts(int c1,int c2){
    //cerr<<"c1 "<<c1<<" c2 "<<endl;
    //taken from CountData.cpp from treemx
    double f = (double) c1 / ( (double) c1 + (double) c2 );
    nInds   += ((double) c1+ (double) c2)/2.0;//number of individuals

    double tmp2 = (double) c2 / ((double) c1+ (double) c2 - 1.0);
    double tmphzy = 2* f * tmp2;

    if((c1+c2) < 2){
	tmphzy = 2*f*(1-f);
    }
    ///cerr<<"c1 "<<c1<<" c2 "<<c2<<" "<<tmphzy<<endl;
    hzy += tmphzy; //2*f*(1-f);
    nsnp++;
}


F3Counter &  F3Counter::operator+=(const F3Counter & other){
    //As f3 is:
    //( (px_1-p1_1)*(px_1-p2_1) +
    //  (px_1-p1_1)*(px_1-p2_1) +
    //  (px_1-p1_1)*(px_1-p2_1) +
    //  (px_1-p1_1)*(px_1-p2_1) )
    //             /
    //            N
    //we should be able to add substract from the f3Sum and counterSites
    this->f3Sum        += other.f3Sum;
    this->counterSites += other.counterSites;
    this->nInds        += other.nInds;
    this->hzy          += other.hzy;
    this->nsnp         += other.nsnp;


    return *this;
}


F3Counter &  F3Counter::operator-=(const F3Counter & other){

    this->f3Sum        -= other.f3Sum;
    this->counterSites -= other.counterSites;
    this->nInds        -= other.nInds;
    this->hzy          -= other.hzy;
    this->nsnp         -= other.nsnp;

    return *this;
}
