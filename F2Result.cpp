/*
 * F2Result
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "F2Result.h"

F2Result::F2Result(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    noDamage.reinitializedCounters();

    // noDamage.reinitializedCounters();
}

F2Result::F2Result(const F2Result & other){
     all             = other.all;
     noCpg           = other.noCpg;
     onlyCpg         = other.onlyCpg;
     transitions     = other.transitions;
     transversions   = other.transversions;
     noDamage        = other.noDamage;
}

F2Result::~F2Result(){

}


string F2Result::printWithBootstrap(list<vector< vector< vector<F2Result> >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int k,unsigned int numberOfBootstraps){
    stringstream toreturn;

    vector<double> allBoot;
    vector<double> nocpgBoot;
    vector<double> onlycpgBoot;
    vector<double> transitionsBoot;
    vector<double> transversionsBoot;
    vector<double> noDamageBoot;

    //all
    for (list< vector< vector< vector<F2Result> >  >   >::iterator it=boostraps.begin(); 
	 it != boostraps.end(); it++){
	allBoot.push_back(              (*it)[i][j][k].all.computeF2()           );
	nocpgBoot.push_back(            (*it)[i][j][k].noCpg.computeF2()         );
	onlycpgBoot.push_back(          (*it)[i][j][k].onlyCpg.computeF2()       );
	transitionsBoot.push_back(      (*it)[i][j][k].transitions.computeF2()   );
	transversionsBoot.push_back(    (*it)[i][j][k].transversions.computeF2() );
	noDamageBoot.push_back(         (*it)[i][j][k].noDamage.computeF2()      );
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
    toreturn<<"all_sites:\t"      <<all           <<"\t"<<allDEV.first           <<"\t"<<allDEV.second           <<"\t"<<all.computeF2()           / allDEV.second <<"\n";
    toreturn<<"nocpg_sites:\t"    <<noCpg         <<"\t"<<noCpgDEV.first         <<"\t"<<noCpgDEV.second         <<"\t"<<noCpg.computeF2()         / noCpgDEV.second <<"\n";
    toreturn<<"onlycpg_sites:\t"  <<onlyCpg       <<"\t"<<onlyCpgDEV.first       <<"\t"<<onlyCpgDEV.second       <<"\t"<<onlyCpg.computeF2()       / onlyCpgDEV.second <<"\n";
    toreturn<<"transitions:\t"    <<transitions   <<"\t"<<transitionsDEV.first   <<"\t"<<transitionsDEV.second   <<"\t"<<transitions.computeF2()   / transitionsDEV.second <<"\n";
    toreturn<<"transversions:\t"  <<transversions <<"\t"<<transversionsDEV.first <<"\t"<<transversionsDEV.second <<"\t"<<transversions.computeF2() / transversionsDEV.second <<"\n";
    toreturn<<"noDamage:\t"       <<noDamage      <<"\t"<<noDamageDEV.first      <<"\t"<<noDamageDEV.second      <<"\t"<<noDamage.computeF2()      / noDamageDEV.second <<"\n";

    return toreturn.str();
}

void F2Result::addIfNotInfinity( vector<double> * target , double val ){
    if(val != std::numeric_limits<double>::infinity()){
	target->push_back(val);
    }
}

string F2Result::printWithJacknife(const vector<const  F2Result * >  * jacknife){
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
	addIfNotInfinity( &allBoot,             jacknife->at(k)->all.computeF2()              );
	addIfNotInfinity( &nocpgBoot,           jacknife->at(k)->noCpg.computeF2()            );
	addIfNotInfinity( &onlycpgBoot,         jacknife->at(k)->onlyCpg.computeF2()          );
	addIfNotInfinity( &transitionsBoot,     jacknife->at(k)->transitions.computeF2()      );
	addIfNotInfinity( &transversionsBoot,   jacknife->at(k)->transversions.computeF2()    );
	addIfNotInfinity( &noDamageBoot,        jacknife->at(k)->noDamage.computeF2()         );


    }


    pair<double,double> allJK           =  computeJackknifeConfIntervals(all.computeF2(),          allBoot          );
    pair<double,double> noCpgJK         =  computeJackknifeConfIntervals(noCpg.computeF2(),        nocpgBoot        );
    pair<double,double> onlyCpgJK       =  computeJackknifeConfIntervals(onlyCpg.computeF2(),      onlycpgBoot      );
    pair<double,double> transitionsJK   =  computeJackknifeConfIntervals(transitions.computeF2(),  transitionsBoot  );
    pair<double,double> transversionsJK =  computeJackknifeConfIntervals(transversions.computeF2(),transversionsBoot);
    pair<double,double> noDamageJK      =  computeJackknifeConfIntervals(noDamage.computeF2(),     transversionsBoot);


    
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

F2Result &  F2Result::operator+=(const F2Result & other){   

    this->all             += other.all;
    this->noCpg           += other.noCpg;
    this->onlyCpg         += other.onlyCpg;
    this->transversions   += other.transversions;
    this->transitions     += other.transitions;
    this->noDamage        += other.noDamage;

    return *this;

}


F2Result &  F2Result::operator-=(const F2Result & other){   
    this->all             -= other.all;
    this->noCpg           -= other.noCpg;
    this->onlyCpg         -= other.onlyCpg;
    this->transversions   -= other.transversions;
    this->transitions     -= other.transitions;
    this->noDamage        -= other.noDamage;

    return *this;

}
