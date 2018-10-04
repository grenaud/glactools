/*
 * F3Result
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "F3Result.h"

F3Result::F3Result(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    noDamage.reinitializedCounters();

    // noDamage.reinitializedCounters();
}

F3Result::F3Result(const F3Result & other){
     all             = other.all;
     noCpg           = other.noCpg;
     onlyCpg         = other.onlyCpg;
     transitions     = other.transitions;
     transversions   = other.transversions;
     noDamage        = other.noDamage;
}

F3Result::~F3Result(){

}


string F3Result::printWithBootstrap(list<vector< vector< vector<F3Result> >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int k,unsigned int numberOfBootstraps){
    stringstream toreturn;

    vector<double> allBoot;
    vector<double> nocpgBoot;
    vector<double> onlycpgBoot;
    vector<double> transitionsBoot;
    vector<double> transversionsBoot;
    vector<double> noDamageBoot;

    //all
    for (list< vector< vector< vector<F3Result> >  >   >::iterator it=boostraps.begin(); 
	 it != boostraps.end(); it++){
	allBoot.push_back(              (*it)[i][j][k].all.computeDST()           );
	nocpgBoot.push_back(            (*it)[i][j][k].noCpg.computeDST()         );
	onlycpgBoot.push_back(          (*it)[i][j][k].onlyCpg.computeDST()       );
	transitionsBoot.push_back(      (*it)[i][j][k].transitions.computeDST()   );
	transversionsBoot.push_back(    (*it)[i][j][k].transversions.computeDST() );
	noDamageBoot.push_back(         (*it)[i][j][k].noDamage.computeDST()      );
    }

    //cout<<vectorToString(allBoot)<<endl;
    pair<double,double> allDEV             = computeMeanSTDDEV(allBoot);
    pair<double,double> noCpgDEV           = computeMeanSTDDEV(nocpgBoot);
    pair<double,double> onlyCpgDEV         = computeMeanSTDDEV(onlycpgBoot);
    pair<double,double> transitionsDEV     = computeMeanSTDDEV(transitionsBoot);
    pair<double,double> transversionsDEV   = computeMeanSTDDEV(transversionsBoot);
    pair<double,double> noDamageDEV        = computeMeanSTDDEV(transversionsBoot);


    //cout<<allBootDEV.first<<"\t"<<allBootDEV.second<<endl;
    toreturn<<":\t"       <<all.headerForCount()  <<"\n";
    toreturn<<"all_sites:\t"      <<all           <<"\t"<<allDEV.first           <<"\t"<<allDEV.second           <<"\t"<<all.computeDST()           / allDEV.second <<"\n";
    toreturn<<"nocpg_sites:\t"    <<noCpg         <<"\t"<<noCpgDEV.first         <<"\t"<<noCpgDEV.second         <<"\t"<<noCpg.computeDST()         / noCpgDEV.second <<"\n";
    toreturn<<"onlycpg_sites:\t"  <<onlyCpg       <<"\t"<<onlyCpgDEV.first       <<"\t"<<onlyCpgDEV.second       <<"\t"<<onlyCpg.computeDST()       / onlyCpgDEV.second <<"\n";
    toreturn<<"transitions:\t"    <<transitions   <<"\t"<<transitionsDEV.first   <<"\t"<<transitionsDEV.second   <<"\t"<<transitions.computeDST()   / transitionsDEV.second <<"\n";
    toreturn<<"transversions:\t"  <<transversions <<"\t"<<transversionsDEV.first <<"\t"<<transversionsDEV.second <<"\t"<<transversions.computeDST() / transversionsDEV.second <<"\n";
    toreturn<<"noDamage:\t"       <<noDamage      <<"\t"<<noDamageDEV.first      <<"\t"<<noDamageDEV.second      <<"\t"<<noDamage.computeDST()      / noDamageDEV.second <<"\n";

    return toreturn.str();
}

void F3Result::addIfNotInfinity( vector<double> * target , double val ){
    if(val != std::numeric_limits<double>::infinity()){
	target->push_back(val);
    }
}

string F3Result::printWithJacknife(const vector<const  F3Result * >  * jacknife){
    stringstream toreturn;

    vector<double> allBoot;
    vector<double> nocpgBoot;
    vector<double> onlycpgBoot;
    vector<double> transitionsBoot;
    vector<double> transversionsBoot;
    vector<double> noDamageBoot;


    // vector<double> counterReferenceB;
    // vector<double> counterBothB;
    
    //all
    //for (vector< AvgCoaResult * > ::iterator it=jacknife->begin(); it != jacknife->end(); it++){
    for (unsigned int k=0;k<jacknife->size(); k++){
	//cout<<"jack "<<jacknife->at(k)->all.avgCoaRefSam()<<endl;
	addIfNotInfinity( &allBoot,             jacknife->at(k)->all.computeDST()              );
	addIfNotInfinity( &nocpgBoot,           jacknife->at(k)->noCpg.computeDST()            );
	addIfNotInfinity( &onlycpgBoot,         jacknife->at(k)->onlyCpg.computeDST()          );
	addIfNotInfinity( &transitionsBoot,     jacknife->at(k)->transitions.computeDST()      );
	addIfNotInfinity( &transversionsBoot,   jacknife->at(k)->transversions.computeDST()    );
	addIfNotInfinity( &noDamageBoot,        jacknife->at(k)->noDamage.computeDST()         );


    }


    pair<double,double> allJK           =  computeJackknifeConfIntervals(all.computeDST(),          allBoot          );
    pair<double,double> noCpgJK         =  computeJackknifeConfIntervals(noCpg.computeDST(),        nocpgBoot        );
    pair<double,double> onlyCpgJK       =  computeJackknifeConfIntervals(onlyCpg.computeDST(),      onlycpgBoot      );
    pair<double,double> transitionsJK   =  computeJackknifeConfIntervals(transitions.computeDST(),  transitionsBoot  );
    pair<double,double> transversionsJK =  computeJackknifeConfIntervals(transversions.computeDST(),transversionsBoot);
    pair<double,double> noDamageJK      =  computeJackknifeConfIntervals(noDamage.computeDST(),     transversionsBoot);


    
    // cout<<"jk size"<<jacknife->size()<<"\t"<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<( (all.avgCoaRefSam()-allDEVRefSam.first)            / allDEVRefSam.second)<<"\t"<<ppp.first<<"\t"<<ppp.second<<"\tlb"<<pppRef.first<<"\tub"<<pppRef.second<<"\t"<<pppBoth.first<<"\t"<<pppBoth.second<<"\t"<<(pppRef.first/pppBoth.second)<<"\t"<<(pppRef.second/pppBoth.first)<<"\t"<<all.counterReference<<"\t"<<(all.counterReference+all.counterCommon)<<endl;
    //exit(1);

    toreturn
	<<all<<"\t"  
	<<allJK.first                <<"\t"<<allJK.second             <<"\t"  

	<<onlyCpg<<"\t"  
	<<onlyCpgJK.first            <<"\t"<<onlyCpgJK.second 	<<"\t"  

	<<noCpg<<"\t"  
	<<noCpgJK.first              <<"\t"<<noCpgJK.second 	        <<"\t"            

	<<transitions<<"\t"  
	<<transitionsJK.first        <<"\t"<<transitionsJK.second 	<<"\t"      

	<<transversions<<"\t"  
	<<transversionsJK.first      <<"\t"<<transversionsJK.second 	<<"\t"    

	<<noDamage<<"\t"  
	<<noDamageJK.first           <<"\t"<<noDamageJK.second       ;
    

    return toreturn.str();

    
}

F3Result &  F3Result::operator+=(const F3Result & other){   

    this->all             += other.all;
    this->noCpg           += other.noCpg;
    this->onlyCpg         += other.onlyCpg;
    this->transversions   += other.transversions;
    this->transitions     += other.transitions;
    this->noDamage        += other.noDamage;

    return *this;

}


F3Result &  F3Result::operator-=(const F3Result & other){   
    this->all             -= other.all;
    this->noCpg           -= other.noCpg;
    this->onlyCpg         -= other.onlyCpg;
    this->transversions   -= other.transversions;
    this->transitions     -= other.transitions;
    this->noDamage        -= other.noDamage;

    return *this;

}
