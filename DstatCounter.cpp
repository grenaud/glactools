/*
 * DstatCounter
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "DstatCounter.h"


DstatCounter::DstatCounter(){
    reinitializedCounters();
}


void DstatCounter::reinitializedCounters(){
    counterAncAnc    =0;
    counterAncDer    =0;
    counterDerAnc    =0;
    counterDerDer    =0;
}


string DstatCounter::headerForCount() const {
    return "AA\tAD\tDA\tDD\t(ADDA-DADA)/(ADDA+ADDA)";
}

vector<double> DstatCounter::performBoostraps() const{
    vector<double> boostraps;
    timeval time;
    gettimeofday(&time, NULL);

    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );

    unsigned int sumADDA = (counterDerAnc+counterAncDer);
    for(unsigned int i=0;i<1000;i++){
	 unsigned int counterAncDerBoot=0;
	 unsigned int counterDerAncBoot=0;
	 for(unsigned int j=0;j<(sumADDA);j++){
	     unsigned int randnum = rand()%sumADDA;
	     if(randnum<counterAncDer)
		 counterAncDerBoot++;
	     else
		 counterDerAncBoot++;
	 }
	 boostraps.push_back( (double(counterAncDerBoot)-double(counterDerAncBoot))/(double(counterAncDerBoot)+double(counterDerAncBoot)) );
    }
    
    return boostraps;
}


double DstatCounter::computeDST() const{
    double dst= (double(counterAncDer)-double(counterDerAnc))/(double(counterAncDer)+double(counterDerAnc));
    return dst;    
}


DstatCounter &  DstatCounter::operator+=(const DstatCounter & other){
    this->counterAncAnc += other.counterAncAnc;
    this->counterAncDer += other.counterAncDer;
    this->counterDerAnc += other.counterDerAnc;
    this->counterDerDer += other.counterDerDer;


    return *this;
}


DstatCounter &  DstatCounter::operator-=(const DstatCounter & other){
    this->counterAncAnc -= other.counterAncAnc;
    this->counterAncDer -= other.counterAncDer;
    this->counterDerAnc -= other.counterDerAnc;
    this->counterDerDer -= other.counterDerDer;


    return *this;
}
