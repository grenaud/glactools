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

vector<double> F3Counter::performBoostraps() const{
    vector<double> boostraps;
    timeval time;
    gettimeofday(&time, NULL);

    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );
    //TODO
    // unsigned int sumADDA = (counterDerAnc+counterAncDer);
    // for(unsigned int i=0;i<1000;i++){
    // 	 unsigned int counterAncDerBoot=0;
    // 	 unsigned int counterDerAncBoot=0;
    // 	 for(unsigned int j=0;j<(sumADDA);j++){
    // 	     unsigned int randnum = rand()%sumADDA;
    // 	     if(randnum<counterAncDer)
    // 		 counterAncDerBoot++;
    // 	     else
    // 		 counterDerAncBoot++;
    // 	 }
    // 	 boostraps.push_back( (double(counterAncDerBoot)-double(counterDerAncBoot))/(double(counterAncDerBoot)+double(counterDerAncBoot)) );
    // }
    
    return boostraps;
}


double F3Counter::computeDST() const{
    //double dst= (double(counterAncDer)-double(counterDerAnc))/(double(counterAncDer)+double(counterDerAnc));
    //TODO
    double dst= (f3Sum)/(double(counterSites));
    return dst;    
}


F3Counter &  F3Counter::operator+=(const F3Counter & other){
    //TODO
    this->f3Sum += other.f3Sum;
    this->counterSites += other.counterSites;


    return *this;
}


F3Counter &  F3Counter::operator-=(const F3Counter & other){
    //TODO
    this->f3Sum        -= other.f3Sum;
    this->counterSites -= other.counterSites;
    // this->counterDerAnc -= other.counterDerAnc;
    // this->counterDerDer -= other.counterDerDer;


    return *this;
}
