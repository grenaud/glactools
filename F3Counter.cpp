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

    f3Sum=0;
    counterSites=0;

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


    return *this;
}


F3Counter &  F3Counter::operator-=(const F3Counter & other){

    this->f3Sum        -= other.f3Sum;
    this->counterSites -= other.counterSites;

    return *this;
}
