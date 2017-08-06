/*
 * AvgCoaResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "AvgCoaResult.h"

AvgCoaResult::AvgCoaResult(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    noDamage.reinitializedCounters();
}

AvgCoaResult::~AvgCoaResult(){

}



string AvgCoaResult::getHeader(){
    //    return "noMut\tcommon\tind1Spec\tind2Spec\tdivInd1\tdivInd1Low\tdivInd1High";
    return 
	all.getHeader("all")+"\t"+
	onlyCpg.getHeader("justCpg")+"\t"+
	noCpg.getHeader("noCpg")+"\t"+
	transitions.getHeader("transi")+"\t"+
	transversions.getHeader("transv")+"\t"+
	noDamage.getHeader("noDam");
}

void AvgCoaResult::addIfNotInfinity( vector<double> * target , double val ){
    if(val != std::numeric_limits<double>::infinity()){
	target->push_back(val);
    }
}

string AvgCoaResult::printWithJacknife(const vector<const  AvgCoaResult * >  * jacknife){
    stringstream toreturn;

    vector<double> allBootRefSam;
    vector<double> nocpgBootRefSam;
    vector<double> onlycpgBootRefSam;
    vector<double> transitionsBootRefSam;
    vector<double> transversionsBootRefSam;
    vector<double> noDamageBootRefSam;

    vector<double> allBootSamRef;
    vector<double> nocpgBootSamRef;
    vector<double> onlycpgBootSamRef;
    vector<double> transitionsBootSamRef;
    vector<double> transversionsBootSamRef;
    vector<double> noDamageBootSamRef;

    // vector<double> counterReferenceB;
    // vector<double> counterBothB;
    
    //all
    //for (vector< AvgCoaResult * > ::iterator it=jacknife->begin(); it != jacknife->end(); it++){
    for (unsigned int k=0;k<jacknife->size(); k++){
	//cout<<"jack "<<jacknife->at(k)->all.avgCoaRefSam()<<endl;
	addIfNotInfinity( &allBootRefSam,             jacknife->at(k)->all.avgCoaRefSam()              );
	addIfNotInfinity( &nocpgBootRefSam,           jacknife->at(k)->noCpg.avgCoaRefSam()            );
	addIfNotInfinity( &onlycpgBootRefSam,         jacknife->at(k)->onlyCpg.avgCoaRefSam()          );
	addIfNotInfinity( &transitionsBootRefSam,     jacknife->at(k)->transitions.avgCoaRefSam()      );
	addIfNotInfinity( &transversionsBootRefSam,   jacknife->at(k)->transversions.avgCoaRefSam()    );
	addIfNotInfinity( &noDamageBootRefSam,        jacknife->at(k)->noDamage.avgCoaRefSam()         );

	addIfNotInfinity( &allBootSamRef,             jacknife->at(k)->all.avgCoaSamRef()           );
	addIfNotInfinity( &nocpgBootSamRef,           jacknife->at(k)->noCpg.avgCoaSamRef()         );
	addIfNotInfinity( &onlycpgBootSamRef,         jacknife->at(k)->onlyCpg.avgCoaSamRef()       );
	addIfNotInfinity( &transitionsBootSamRef,     jacknife->at(k)->transitions.avgCoaSamRef()   );
	addIfNotInfinity( &transversionsBootSamRef,   jacknife->at(k)->transversions.avgCoaSamRef() );
	addIfNotInfinity( &noDamageBootSamRef,        jacknife->at(k)->noDamage.avgCoaSamRef()      );

	// counterReferenceB.push_back(jacknife->at(k)->all.counterReference);
	// counterBothB.push_back(jacknife->at(k)->all.counterReference+jacknife->at(k)->all.counterCommon);

    }

    // cout<<"jk size"<<jacknife->size()<<endl;

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

    pair<double,double> allDEVRefSamJK           =  computeJackknifeConfIntervals(all.avgCoaRefSam(),allBootRefSam);
    pair<double,double> noCpgDEVRefSamJK         =  computeJackknifeConfIntervals(noCpg.avgCoaRefSam(),nocpgBootRefSam);
    pair<double,double> onlyCpgDEVRefSamJK       =  computeJackknifeConfIntervals(onlyCpg.avgCoaRefSam(),onlycpgBootRefSam);
    pair<double,double> transitionsDEVRefSamJK   =  computeJackknifeConfIntervals(transitions.avgCoaRefSam(),transitionsBootRefSam);
    pair<double,double> transversionsDEVRefSamJK =  computeJackknifeConfIntervals(transversions.avgCoaRefSam(),transversionsBootRefSam);
    pair<double,double> noDamageDEVRefSamJK      =  computeJackknifeConfIntervals(noDamage.avgCoaRefSam(),transversionsBootRefSam);

    pair<double,double> allDEVSamRefJK           =  computeJackknifeConfIntervals(all.avgCoaRefSam(),allBootSamRef);
    pair<double,double> noCpgDEVSamRefJK         =  computeJackknifeConfIntervals(noCpg.avgCoaRefSam(),nocpgBootSamRef);
    pair<double,double> onlyCpgDEVSamRefJK       =  computeJackknifeConfIntervals(onlyCpg.avgCoaRefSam(),onlycpgBootSamRef);
    pair<double,double> transitionsDEVSamRefJK   =  computeJackknifeConfIntervals(transitions.avgCoaRefSam(),transitionsBootSamRef);
    pair<double,double> transversionsDEVSamRefJK =  computeJackknifeConfIntervals(transversions.avgCoaRefSam(),transversionsBootSamRef);
    pair<double,double> noDamageDEVSamRefJK      =  computeJackknifeConfIntervals(noDamage.avgCoaRefSam(),transversionsBootSamRef);


    // pair<double,double> ppp     = computeJackknifeConfIntervals(all.avgCoaRefSam(),                        allBootRefSam);
    // pair<double,double> pppRef  = computeJackknifeConfIntervals(all.counterReference,                   counterReferenceB);
    // pair<double,double> pppBoth = computeJackknifeConfIntervals(all.counterReference+all.counterCommon, counterBothB);

    
    // cout<<"jk size"<<jacknife->size()<<"\t"<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<( (all.avgCoaRefSam()-allDEVRefSam.first)            / allDEVRefSam.second)<<"\t"<<ppp.first<<"\t"<<ppp.second<<"\tlb"<<pppRef.first<<"\tub"<<pppRef.second<<"\t"<<pppBoth.first<<"\t"<<pppBoth.second<<"\t"<<(pppRef.first/pppBoth.second)<<"\t"<<(pppRef.second/pppBoth.first)<<"\t"<<all.counterReference<<"\t"<<(all.counterReference+all.counterCommon)<<endl;
    //exit(1);

    toreturn
	<<all<<"\t"  
	<<allDEVRefSamJK.first                <<"\t"<<allDEVRefSamJK.second             <<"\t"  
	<<allDEVSamRefJK.first                <<"\t"<<allDEVSamRefJK.second             <<"\t"  

	<<onlyCpg<<"\t"  
	<<onlyCpgDEVRefSamJK.first            <<"\t"<<onlyCpgDEVRefSamJK.second 	<<"\t"  
	<<onlyCpgDEVSamRefJK.first            <<"\t"<<onlyCpgDEVSamRefJK.second 	<<"\t"          

	<<noCpg<<"\t"  
	<<noCpgDEVRefSamJK.first              <<"\t"<<noCpgDEVRefSamJK.second 	        <<"\t"            
	<<noCpgDEVSamRefJK.first              <<"\t"<<noCpgDEVSamRefJK.second     	<<"\t"        

	<<transitions<<"\t"  
	<<transitionsDEVRefSamJK.first        <<"\t"<<transitionsDEVRefSamJK.second 	<<"\t"      
	<<transitionsDEVSamRefJK.first        <<"\t"<<transitionsDEVSamRefJK.second 	<<"\t"      

	<<transversions<<"\t"  
	<<transversionsDEVRefSamJK.first      <<"\t"<<transversionsDEVRefSamJK.second 	<<"\t"    
	<<transversionsDEVSamRefJK.first      <<"\t"<<transversionsDEVSamRefJK.second 	<<"\t"    

	<<noDamage<<"\t"  
	<<noDamageDEVRefSamJK.first           <<"\t"<<noDamageDEVRefSamJK.second   	<<"\t"       
	<<noDamageDEVSamRefJK.first           <<"\t"<<noDamageDEVSamRefJK.second       ;
    

    return toreturn.str();

    
}

string AvgCoaResult::printWithBootstrap(list< vector< vector< AvgCoaResult >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int numberOfBootstraps){
    stringstream toreturn;

    vector<double> allBootRefSam;
    vector<double> nocpgBootRefSam;
    vector<double> onlycpgBootRefSam;
    vector<double> transitionsBootRefSam;
    vector<double> transversionsBootRefSam;
    vector<double> noDamageBootRefSam;

    vector<double> allBootSamRef;
    vector<double> nocpgBootSamRef;
    vector<double> onlycpgBootSamRef;
    vector<double> transitionsBootSamRef;
    vector<double> transversionsBootSamRef;
    vector<double> noDamageBootSamRef;

    //all
    for (list< vector< vector< AvgCoaResult >  >   >::iterator it=boostraps.begin(); 
	 it != boostraps.end(); it++){


	addIfNotInfinity( &allBootRefSam,             (*it)[i][j].all.avgCoaRefSam()              );
	addIfNotInfinity( &nocpgBootRefSam,           (*it)[i][j].noCpg.avgCoaRefSam()            );
	addIfNotInfinity( &onlycpgBootRefSam,         (*it)[i][j].onlyCpg.avgCoaRefSam()          );
	addIfNotInfinity( &transitionsBootRefSam,     (*it)[i][j].transitions.avgCoaRefSam()      );
	addIfNotInfinity( &transversionsBootRefSam,   (*it)[i][j].transversions.avgCoaRefSam()    );
	addIfNotInfinity( &noDamageBootRefSam,        (*it)[i][j].noDamage.avgCoaRefSam()         );

	addIfNotInfinity( &allBootSamRef,             (*it)[i][j].all.avgCoaSamRef()           );
	addIfNotInfinity( &nocpgBootSamRef,           (*it)[i][j].noCpg.avgCoaSamRef()         );
	addIfNotInfinity( &onlycpgBootSamRef,         (*it)[i][j].onlyCpg.avgCoaSamRef()       );
	addIfNotInfinity( &transitionsBootSamRef,     (*it)[i][j].transitions.avgCoaSamRef()   );
	addIfNotInfinity( &transversionsBootSamRef,   (*it)[i][j].transversions.avgCoaSamRef() );
	addIfNotInfinity( &noDamageBootSamRef,        (*it)[i][j].noDamage.avgCoaSamRef()      );

					    
    }

    pair<double,double> allDEVRefSam           =  computeMeanSTDDEV(allBootRefSam);
    pair<double,double> noCpgDEVRefSam         =  computeMeanSTDDEV(nocpgBootRefSam);
    pair<double,double> onlyCpgDEVRefSam       =  computeMeanSTDDEV(onlycpgBootRefSam);
    pair<double,double> transitionsDEVRefSam   =  computeMeanSTDDEV(transitionsBootRefSam);
    pair<double,double> transversionsDEVRefSam =  computeMeanSTDDEV(transversionsBootRefSam);
    pair<double,double> noDamageDEVRefSam      =  computeMeanSTDDEV(transversionsBootRefSam);

    pair<double,double> allDEVSamRef           =  computeMeanSTDDEV(allBootSamRef);
    pair<double,double> noCpgDEVSamRef         =  computeMeanSTDDEV(nocpgBootSamRef);
    pair<double,double> onlyCpgDEVSamRef       =  computeMeanSTDDEV(onlycpgBootSamRef);
    pair<double,double> transitionsDEVSamRef   =  computeMeanSTDDEV(transitionsBootSamRef);
    pair<double,double> transversionsDEVSamRef =  computeMeanSTDDEV(transversionsBootSamRef);
    pair<double,double> noDamageDEVSamRef      =  computeMeanSTDDEV(transversionsBootSamRef);

    toreturn
	<<all<<"\t"  
	<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<(all.avgCoaRefSam()            / allDEVRefSam.second) <<"\t"
	<<allDEVSamRef.first                <<"\t"<<allDEVSamRef.second            <<"\t"<<(all.avgCoaSamRef()            / allDEVSamRef.second) <<"\t"

	<<onlyCpg<<"\t"  
	<<onlyCpgDEVRefSam.first            <<"\t"<<onlyCpgDEVRefSam.second        <<"\t"<<(onlyCpg.avgCoaRefSam()        / onlyCpgDEVRefSam.second) <<"\t"
	<<onlyCpgDEVSamRef.first            <<"\t"<<onlyCpgDEVSamRef.second        <<"\t"<<(onlyCpg.avgCoaSamRef()        / onlyCpgDEVSamRef.second) <<"\t"

	<<noCpg<<"\t"  
	<<noCpgDEVRefSam.first              <<"\t"<<noCpgDEVRefSam.second          <<"\t"<<(noCpg.avgCoaRefSam()          / noCpgDEVRefSam.second) <<"\t"
	<<noCpgDEVSamRef.first              <<"\t"<<noCpgDEVSamRef.second          <<"\t"<<(noCpg.avgCoaSamRef()          / noCpgDEVSamRef.second) <<"\t"

	<<transitions<<"\t"  
	<<transitionsDEVRefSam.first        <<"\t"<<transitionsDEVRefSam.second    <<"\t"<<(transitions.avgCoaRefSam()    / transitionsDEVRefSam.second) <<"\t"
	<<transitionsDEVSamRef.first        <<"\t"<<transitionsDEVSamRef.second    <<"\t"<<(transitions.avgCoaSamRef()    / transitionsDEVSamRef.second) <<"\t"

	<<transversions<<"\t"  
	<<transversionsDEVRefSam.first      <<"\t"<<transversionsDEVRefSam.second  <<"\t"<<(transversions.avgCoaRefSam()  / transversionsDEVRefSam.second) <<"\t"
	<<transversionsDEVSamRef.first      <<"\t"<<transversionsDEVSamRef.second  <<"\t"<<(transversions.avgCoaSamRef()  / transversionsDEVSamRef.second) <<"\t"

	<<noDamage<<"\t"  
	<<noDamageDEVRefSam.first           <<"\t"<<noDamageDEVRefSam.second       <<"\t"<<(noDamage.avgCoaRefSam()       / noDamageDEVRefSam.second) <<"\t"
	<<noDamageDEVSamRef.first           <<"\t"<<noDamageDEVSamRef.second       <<"\t"<<(noDamage.avgCoaSamRef()       / noDamageDEVSamRef.second) ;
    

    return toreturn.str();
    
}
