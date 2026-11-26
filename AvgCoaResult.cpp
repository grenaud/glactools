/*
 * AvgCoaResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "AvgCoaResult.h"
#include "SafeJackknife.h"

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

void AvgCoaResult::addReplicate( vector<double> * target , double val ){
    target->push_back(val);
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
	addReplicate( &allBootRefSam,             jacknife->at(k)->all.avgCoaRefSam()              );
	addReplicate( &nocpgBootRefSam,           jacknife->at(k)->noCpg.avgCoaRefSam()            );
	addReplicate( &onlycpgBootRefSam,         jacknife->at(k)->onlyCpg.avgCoaRefSam()          );
	addReplicate( &transitionsBootRefSam,     jacknife->at(k)->transitions.avgCoaRefSam()      );
	addReplicate( &transversionsBootRefSam,   jacknife->at(k)->transversions.avgCoaRefSam()    );
	addReplicate( &noDamageBootRefSam,        jacknife->at(k)->noDamage.avgCoaRefSam()         );

	addReplicate( &allBootSamRef,             jacknife->at(k)->all.avgCoaSamRef()           );
	addReplicate( &nocpgBootSamRef,           jacknife->at(k)->noCpg.avgCoaSamRef()         );
	addReplicate( &onlycpgBootSamRef,         jacknife->at(k)->onlyCpg.avgCoaSamRef()       );
	addReplicate( &transitionsBootSamRef,     jacknife->at(k)->transitions.avgCoaSamRef()   );
	addReplicate( &transversionsBootSamRef,   jacknife->at(k)->transversions.avgCoaSamRef() );
	addReplicate( &noDamageBootSamRef,        jacknife->at(k)->noDamage.avgCoaSamRef()      );

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

    pair<double,double> allDEVRefSamJK           =  jackknifeOrNan(all.avgCoaRefSam(),allBootRefSam);
    pair<double,double> noCpgDEVRefSamJK         =  jackknifeOrNan(noCpg.avgCoaRefSam(),nocpgBootRefSam);
    pair<double,double> onlyCpgDEVRefSamJK       =  jackknifeOrNan(onlyCpg.avgCoaRefSam(),onlycpgBootRefSam);
    pair<double,double> transitionsDEVRefSamJK   =  jackknifeOrNan(transitions.avgCoaRefSam(),transitionsBootRefSam);
    pair<double,double> transversionsDEVRefSamJK =  jackknifeOrNan(transversions.avgCoaRefSam(),transversionsBootRefSam);
    pair<double,double> noDamageDEVRefSamJK      =  jackknifeOrNan(noDamage.avgCoaRefSam(),noDamageBootRefSam);

    pair<double,double> allDEVSamRefJK           =  jackknifeOrNan(all.avgCoaSamRef(),allBootSamRef);
    pair<double,double> noCpgDEVSamRefJK         =  jackknifeOrNan(noCpg.avgCoaSamRef(),nocpgBootSamRef);
    pair<double,double> onlyCpgDEVSamRefJK       =  jackknifeOrNan(onlyCpg.avgCoaSamRef(),onlycpgBootSamRef);
    pair<double,double> transitionsDEVSamRefJK   =  jackknifeOrNan(transitions.avgCoaSamRef(),transitionsBootSamRef);
    pair<double,double> transversionsDEVSamRefJK =  jackknifeOrNan(transversions.avgCoaSamRef(),transversionsBootSamRef);
    pair<double,double> noDamageDEVSamRefJK      =  jackknifeOrNan(noDamage.avgCoaSamRef(),noDamageBootSamRef);


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


	addReplicate( &allBootRefSam,             (*it)[i][j].all.avgCoaRefSam()              );
	addReplicate( &nocpgBootRefSam,           (*it)[i][j].noCpg.avgCoaRefSam()            );
	addReplicate( &onlycpgBootRefSam,         (*it)[i][j].onlyCpg.avgCoaRefSam()          );
	addReplicate( &transitionsBootRefSam,     (*it)[i][j].transitions.avgCoaRefSam()      );
	addReplicate( &transversionsBootRefSam,   (*it)[i][j].transversions.avgCoaRefSam()    );
	addReplicate( &noDamageBootRefSam,        (*it)[i][j].noDamage.avgCoaRefSam()         );

	addReplicate( &allBootSamRef,             (*it)[i][j].all.avgCoaSamRef()           );
	addReplicate( &nocpgBootSamRef,           (*it)[i][j].noCpg.avgCoaSamRef()         );
	addReplicate( &onlycpgBootSamRef,         (*it)[i][j].onlyCpg.avgCoaSamRef()       );
	addReplicate( &transitionsBootSamRef,     (*it)[i][j].transitions.avgCoaSamRef()   );
	addReplicate( &transversionsBootSamRef,   (*it)[i][j].transversions.avgCoaSamRef() );
	addReplicate( &noDamageBootSamRef,        (*it)[i][j].noDamage.avgCoaSamRef()      );

					    
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
