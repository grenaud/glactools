/*
 * DistResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "DistResult.h"

DistResult::DistResult(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    noDamage.reinitializedCounters();
}

DistResult::~DistResult(){

}



string DistResult::getHeader(){
    //    return "noMut\tcommon\tind1Spec\tind2Spec\tdivInd1\tdivInd1Low\tdivInd1High";
    return 
	all.getHeader("all")+"\t"+
	onlyCpg.getHeader("justCpg")+"\t"+
	noCpg.getHeader("noCpg")+"\t"+
	transitions.getHeader("transi")+"\t"+
	transversions.getHeader("transv")+"\t"+
	noDamage.getHeader("noDam");
}

void DistResult::addIfNotInfinity( vector<double> * target , double val ){
    if(val != std::numeric_limits<double>::infinity()){
	target->push_back(val);
    }
}

string DistResult::printWithJacknife(const vector<const  DistResult * >  * jacknife,const string & dnaDistMode){
    stringstream toreturn;
    //cerr<<"printWithJacknife()"<<endl;
    vector<double> allBootDist;
    vector<double> nocpgBootDist;
    vector<double> onlycpgBootDist;
    vector<double> transitionsBootDist;
    vector<double> transversionsBootDist;
    vector<double> noDamageBootDist;

    // vector<double> allBootSamRef;
    // vector<double> nocpgBootSamRef;
    // vector<double> onlycpgBootSamRef;
    // vector<double> transitionsBootSamRef;
    // vector<double> transversionsBootSamRef;
    // vector<double> noDamageBootSamRef;

    // // vector<double> counterReferenceB;
    // // vector<double> counterBothB;
    
    // //all
    // //for (vector< DistResult * > ::iterator it=jacknife->begin(); it != jacknife->end(); it++){
    for (unsigned int k=0;k<jacknife->size(); k++){
    	//cout<<"jack "<<jacknife->at(k)->all.distRefSam()<<endl;
    	addIfNotInfinity( &allBootDist,             jacknife->at(k)->all.dist(          dnaDistMode)    );
    	addIfNotInfinity( &nocpgBootDist,           jacknife->at(k)->noCpg.dist(        dnaDistMode)    );
    	addIfNotInfinity( &onlycpgBootDist,         jacknife->at(k)->onlyCpg.dist(      dnaDistMode)    );
    	addIfNotInfinity( &transitionsBootDist,     jacknife->at(k)->transitions.dist(  dnaDistMode)    );
    	addIfNotInfinity( &transversionsBootDist,   jacknife->at(k)->transversions.dist(dnaDistMode)    );
    	addIfNotInfinity( &noDamageBootDist,        jacknife->at(k)->noDamage.dist(     dnaDistMode)    );

    }
    // }

    // // cout<<"jk size"<<jacknife->size()<<endl;

    // pair<double,double> allDEVDist           =  computeMeanSTDDEV(allBootDist);
    // pair<double,double> noCpgDEVDist         =  computeMeanSTDDEV(nocpgBootDist);
    // pair<double,double> onlyCpgDEVDist       =  computeMeanSTDDEV(onlycpgBootDist);
    // pair<double,double> transitionsDEVDist   =  computeMeanSTDDEV(transitionsBootDist);
    // pair<double,double> transversionsDEVDist =  computeMeanSTDDEV(transversionsBootDist);
    // pair<double,double> noDamageDEVDist      =  computeMeanSTDDEV(transversionsBootDist);

    // // pair<double,double> allDEVSamRef           =  computeMeanSTDDEV(allBootSamRef);
    // // pair<double,double> noCpgDEVSamRef         =  computeMeanSTDDEV(nocpgBootSamRef);
    // // pair<double,double> onlyCpgDEVSamRef       =  computeMeanSTDDEV(onlycpgBootSamRef);
    // // pair<double,double> transitionsDEVSamRef   =  computeMeanSTDDEV(transitionsBootSamRef);
    // // pair<double,double> transversionsDEVSamRef =  computeMeanSTDDEV(transversionsBootSamRef);
    // // pair<double,double> noDamageDEVSamRef      =  computeMeanSTDDEV(transversionsBootSamRef);

    pair<double,double> allDEVDistJK           =  computeJackknifeConfIntervals(all.dist(dnaDistMode),          allBootDist);
    pair<double,double> noCpgDEVDistJK         =  computeJackknifeConfIntervals(noCpg.dist(dnaDistMode),        nocpgBootDist);
    pair<double,double> onlyCpgDEVDistJK       =  computeJackknifeConfIntervals(onlyCpg.dist(dnaDistMode),      onlycpgBootDist);
    pair<double,double> transitionsDEVDistJK   =  computeJackknifeConfIntervals(transitions.dist(dnaDistMode),  transitionsBootDist);
    pair<double,double> transversionsDEVDistJK =  computeJackknifeConfIntervals(transversions.dist(dnaDistMode),transversionsBootDist);
    pair<double,double> noDamageDEVDistJK      =  computeJackknifeConfIntervals(noDamage.dist(dnaDistMode),     transversionsBootDist);

    // pair<double,double> allDEVSamRefJK           =  computeJackknifeConfIntervals(all.distRefSam(),allBootSamRef);
    // pair<double,double> noCpgDEVSamRefJK         =  computeJackknifeConfIntervals(noCpg.distRefSam(),nocpgBootSamRef);
    // pair<double,double> onlyCpgDEVSamRefJK       =  computeJackknifeConfIntervals(onlyCpg.distRefSam(),onlycpgBootSamRef);
    // pair<double,double> transitionsDEVSamRefJK   =  computeJackknifeConfIntervals(transitions.distRefSam(),transitionsBootSamRef);
    // pair<double,double> transversionsDEVSamRefJK =  computeJackknifeConfIntervals(transversions.distRefSam(),transversionsBootSamRef);
    // pair<double,double> noDamageDEVSamRefJK      =  computeJackknifeConfIntervals(noDamage.distRefSam(),transversionsBootSamRef);


    // // pair<double,double> ppp     = computeJackknifeConfIntervals(all.distRefSam(),                        allBootRefSam);
    // // pair<double,double> pppRef  = computeJackknifeConfIntervals(all.counterReference,                   counterReferenceB);
    // // pair<double,double> pppBoth = computeJackknifeConfIntervals(all.counterReference+all.counterCommon, counterBothB);

    
    // // cout<<"jk size"<<jacknife->size()<<"\t"<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<( (all.distRefSam()-allDEVRefSam.first)            / allDEVRefSam.second)<<"\t"<<ppp.first<<"\t"<<ppp.second<<"\tlb"<<pppRef.first<<"\tub"<<pppRef.second<<"\t"<<pppBoth.first<<"\t"<<pppBoth.second<<"\t"<<(pppRef.first/pppBoth.second)<<"\t"<<(pppRef.second/pppBoth.first)<<"\t"<<all.counterReference<<"\t"<<(all.counterReference+all.counterCommon)<<endl;
    // //exit(1);

    toreturn
    	<<all<<"\t"  
    	<<allDEVDistJK.first                <<"\t"<<allDEVDistJK.second             <<"\t"  

    	<<onlyCpg<<"\t"  
    	<<onlyCpgDEVDistJK.first            <<"\t"<<onlyCpgDEVDistJK.second 	<<"\t"  

    	<<noCpg<<"\t"  
    	<<noCpgDEVDistJK.first              <<"\t"<<noCpgDEVDistJK.second 	        <<"\t"            

    	<<transitions<<"\t"  
    	<<transitionsDEVDistJK.first        <<"\t"<<transitionsDEVDistJK.second 	<<"\t"      

    	<<transversions<<"\t"  
    	<<transversionsDEVDistJK.first      <<"\t"<<transversionsDEVDistJK.second 	<<"\t"    

    	<<noDamage<<"\t"  
    	<<noDamageDEVDistJK.first           <<"\t"<<noDamageDEVDistJK.second   	;
    return toreturn.str();

    
}

string DistResult::printWithBootstrap(list< vector< vector< DistResult >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int numberOfBootstraps){
    stringstream toreturn;

    // vector<double> allBootRefSam;
    // vector<double> nocpgBootRefSam;
    // vector<double> onlycpgBootRefSam;
    // vector<double> transitionsBootRefSam;
    // vector<double> transversionsBootRefSam;
    // vector<double> noDamageBootRefSam;

    // vector<double> allBootSamRef;
    // vector<double> nocpgBootSamRef;
    // vector<double> onlycpgBootSamRef;
    // vector<double> transitionsBootSamRef;
    // vector<double> transversionsBootSamRef;
    // vector<double> noDamageBootSamRef;

    // //all
    // for (list< vector< vector< DistResult >  >   >::iterator it=boostraps.begin(); 
    // 	 it != boostraps.end(); it++){


    // 	addIfNotInfinity( &allBootRefSam,             (*it)[i][j].all.distRefSam()              );
    // 	addIfNotInfinity( &nocpgBootRefSam,           (*it)[i][j].noCpg.distRefSam()            );
    // 	addIfNotInfinity( &onlycpgBootRefSam,         (*it)[i][j].onlyCpg.distRefSam()          );
    // 	addIfNotInfinity( &transitionsBootRefSam,     (*it)[i][j].transitions.distRefSam()      );
    // 	addIfNotInfinity( &transversionsBootRefSam,   (*it)[i][j].transversions.distRefSam()    );
    // 	addIfNotInfinity( &noDamageBootRefSam,        (*it)[i][j].noDamage.distRefSam()         );

    // 	addIfNotInfinity( &allBootSamRef,             (*it)[i][j].all.distSamRef()           );
    // 	addIfNotInfinity( &nocpgBootSamRef,           (*it)[i][j].noCpg.distSamRef()         );
    // 	addIfNotInfinity( &onlycpgBootSamRef,         (*it)[i][j].onlyCpg.distSamRef()       );
    // 	addIfNotInfinity( &transitionsBootSamRef,     (*it)[i][j].transitions.distSamRef()   );
    // 	addIfNotInfinity( &transversionsBootSamRef,   (*it)[i][j].transversions.distSamRef() );
    // 	addIfNotInfinity( &noDamageBootSamRef,        (*it)[i][j].noDamage.distSamRef()      );

					    
    // }

    // pair<double,double> allDEVRefSam           =  computeMeanSTDDEV(allBootRefSam);
    // pair<double,double> noCpgDEVRefSam         =  computeMeanSTDDEV(nocpgBootRefSam);
    // pair<double,double> onlyCpgDEVRefSam       =  computeMeanSTDDEV(onlycpgBootRefSam);
    // pair<double,double> transitionsDEVRefSam   =  computeMeanSTDDEV(transitionsBootRefSam);
    // pair<double,double> transversionsDEVRefSam =  computeMeanSTDDEV(transversionsBootRefSam);
    // pair<double,double> noDamageDEVRefSam      =  computeMeanSTDDEV(transversionsBootRefSam);

    // pair<double,double> allDEVSamRef           =  computeMeanSTDDEV(allBootSamRef);
    // pair<double,double> noCpgDEVSamRef         =  computeMeanSTDDEV(nocpgBootSamRef);
    // pair<double,double> onlyCpgDEVSamRef       =  computeMeanSTDDEV(onlycpgBootSamRef);
    // pair<double,double> transitionsDEVSamRef   =  computeMeanSTDDEV(transitionsBootSamRef);
    // pair<double,double> transversionsDEVSamRef =  computeMeanSTDDEV(transversionsBootSamRef);
    // pair<double,double> noDamageDEVSamRef      =  computeMeanSTDDEV(transversionsBootSamRef);

    // toreturn
    // 	<<all<<"\t"  
    // 	<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<(all.distRefSam()            / allDEVRefSam.second) <<"\t"
    // 	<<allDEVSamRef.first                <<"\t"<<allDEVSamRef.second            <<"\t"<<(all.distSamRef()            / allDEVSamRef.second) <<"\t"

    // 	<<onlyCpg<<"\t"  
    // 	<<onlyCpgDEVRefSam.first            <<"\t"<<onlyCpgDEVRefSam.second        <<"\t"<<(onlyCpg.distRefSam()        / onlyCpgDEVRefSam.second) <<"\t"
    // 	<<onlyCpgDEVSamRef.first            <<"\t"<<onlyCpgDEVSamRef.second        <<"\t"<<(onlyCpg.distSamRef()        / onlyCpgDEVSamRef.second) <<"\t"

    // 	<<noCpg<<"\t"  
    // 	<<noCpgDEVRefSam.first              <<"\t"<<noCpgDEVRefSam.second          <<"\t"<<(noCpg.distRefSam()          / noCpgDEVRefSam.second) <<"\t"
    // 	<<noCpgDEVSamRef.first              <<"\t"<<noCpgDEVSamRef.second          <<"\t"<<(noCpg.distSamRef()          / noCpgDEVSamRef.second) <<"\t"

    // 	<<transitions<<"\t"  
    // 	<<transitionsDEVRefSam.first        <<"\t"<<transitionsDEVRefSam.second    <<"\t"<<(transitions.distRefSam()    / transitionsDEVRefSam.second) <<"\t"
    // 	<<transitionsDEVSamRef.first        <<"\t"<<transitionsDEVSamRef.second    <<"\t"<<(transitions.distSamRef()    / transitionsDEVSamRef.second) <<"\t"

    // 	<<transversions<<"\t"  
    // 	<<transversionsDEVRefSam.first      <<"\t"<<transversionsDEVRefSam.second  <<"\t"<<(transversions.distRefSam()  / transversionsDEVRefSam.second) <<"\t"
    // 	<<transversionsDEVSamRef.first      <<"\t"<<transversionsDEVSamRef.second  <<"\t"<<(transversions.distSamRef()  / transversionsDEVSamRef.second) <<"\t"

    // 	<<noDamage<<"\t"  
    // 	<<noDamageDEVRefSam.first           <<"\t"<<noDamageDEVRefSam.second       <<"\t"<<(noDamage.distRefSam()       / noDamageDEVRefSam.second) <<"\t"
    // 	<<noDamageDEVSamRef.first           <<"\t"<<noDamageDEVSamRef.second       <<"\t"<<(noDamage.distSamRef()       / noDamageDEVSamRef.second) ;
    

    return toreturn.str();
    
}


string DistResult::printAll() const{
    stringstream toreturn;
    toreturn<<all.printWithLabel();
    //toreturn<<all;
    return toreturn.str();
}
