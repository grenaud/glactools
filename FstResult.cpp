/*
 * FstResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FstResult.h"
#include "SafeJackknife.h"

FstResult::FstResult(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    noDamage.reinitializedCounters();
}

FstResult::~FstResult(){

}



string FstResult::getHeader(){
    //    return "noMut\tcommon\tind1Spec\tind2Spec\tdivInd1\tdivInd1Low\tdivInd1High";
    return 
	all.getHeader("all")+"\t"+
	onlyCpg.getHeader("justCpg")+"\t"+
	noCpg.getHeader("noCpg")+"\t"+
	transitions.getHeader("transi")+"\t"+
	transversions.getHeader("transv")+"\t"+
	noDamage.getHeader("noDam");
}

void FstResult::addReplicate( vector<double> * target , double val ){
    target->push_back(val);
}

string FstResult::printWithJacknife(const vector<const  FstResult * >  * jacknife){
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
    //for (vector< FstResult * > ::iterator it=jacknife->begin(); it != jacknife->end(); it++){
    for (unsigned int k=0;k<jacknife->size(); k++){
	//cout<<"jack "<<jacknife->at(k)->all.fstRefSam()<<endl;
	addReplicate( &allBootRefSam,             jacknife->at(k)->all.fstRefSam()              );
	addReplicate( &nocpgBootRefSam,           jacknife->at(k)->noCpg.fstRefSam()            );
	addReplicate( &onlycpgBootRefSam,         jacknife->at(k)->onlyCpg.fstRefSam()          );
	addReplicate( &transitionsBootRefSam,     jacknife->at(k)->transitions.fstRefSam()      );
	addReplicate( &transversionsBootRefSam,   jacknife->at(k)->transversions.fstRefSam()    );
	addReplicate( &noDamageBootRefSam,        jacknife->at(k)->noDamage.fstRefSam()         );

	addReplicate( &allBootSamRef,             jacknife->at(k)->all.fstSamRef()           );
	addReplicate( &nocpgBootSamRef,           jacknife->at(k)->noCpg.fstSamRef()         );
	addReplicate( &onlycpgBootSamRef,         jacknife->at(k)->onlyCpg.fstSamRef()       );
	addReplicate( &transitionsBootSamRef,     jacknife->at(k)->transitions.fstSamRef()   );
	addReplicate( &transversionsBootSamRef,   jacknife->at(k)->transversions.fstSamRef() );
	addReplicate( &noDamageBootSamRef,        jacknife->at(k)->noDamage.fstSamRef()      );

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

    pair<double,double> allDEVRefSamJK           =  jackknifeOrNan(all.fstRefSam(),allBootRefSam);
    pair<double,double> noCpgDEVRefSamJK         =  jackknifeOrNan(noCpg.fstRefSam(),nocpgBootRefSam);
    pair<double,double> onlyCpgDEVRefSamJK       =  jackknifeOrNan(onlyCpg.fstRefSam(),onlycpgBootRefSam);
    pair<double,double> transitionsDEVRefSamJK   =  jackknifeOrNan(transitions.fstRefSam(),transitionsBootRefSam);
    pair<double,double> transversionsDEVRefSamJK =  jackknifeOrNan(transversions.fstRefSam(),transversionsBootRefSam);
    pair<double,double> noDamageDEVRefSamJK      =  jackknifeOrNan(noDamage.fstRefSam(),transversionsBootRefSam);

    pair<double,double> allDEVSamRefJK           =  jackknifeOrNan(all.fstRefSam(),allBootSamRef);
    pair<double,double> noCpgDEVSamRefJK         =  jackknifeOrNan(noCpg.fstRefSam(),nocpgBootSamRef);
    pair<double,double> onlyCpgDEVSamRefJK       =  jackknifeOrNan(onlyCpg.fstRefSam(),onlycpgBootSamRef);
    pair<double,double> transitionsDEVSamRefJK   =  jackknifeOrNan(transitions.fstRefSam(),transitionsBootSamRef);
    pair<double,double> transversionsDEVSamRefJK =  jackknifeOrNan(transversions.fstRefSam(),transversionsBootSamRef);
    pair<double,double> noDamageDEVSamRefJK      =  jackknifeOrNan(noDamage.fstRefSam(),transversionsBootSamRef);


    // pair<double,double> ppp     = computeJackknifeConfIntervals(all.fstRefSam(),                        allBootRefSam);
    // pair<double,double> pppRef  = computeJackknifeConfIntervals(all.counterReference,                   counterReferenceB);
    // pair<double,double> pppBoth = computeJackknifeConfIntervals(all.counterReference+all.counterCommon, counterBothB);

    
    // cout<<"jk size"<<jacknife->size()<<"\t"<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<( (all.fstRefSam()-allDEVRefSam.first)            / allDEVRefSam.second)<<"\t"<<ppp.first<<"\t"<<ppp.second<<"\tlb"<<pppRef.first<<"\tub"<<pppRef.second<<"\t"<<pppBoth.first<<"\t"<<pppBoth.second<<"\t"<<(pppRef.first/pppBoth.second)<<"\t"<<(pppRef.second/pppBoth.first)<<"\t"<<all.counterReference<<"\t"<<(all.counterReference+all.counterCommon)<<endl;
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

string FstResult::printWithBootstrap(list< vector< vector< FstResult >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int numberOfBootstraps){
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
    for (list< vector< vector< FstResult >  >   >::iterator it=boostraps.begin(); 
	 it != boostraps.end(); it++){


	addReplicate( &allBootRefSam,             (*it)[i][j].all.fstRefSam()              );
	addReplicate( &nocpgBootRefSam,           (*it)[i][j].noCpg.fstRefSam()            );
	addReplicate( &onlycpgBootRefSam,         (*it)[i][j].onlyCpg.fstRefSam()          );
	addReplicate( &transitionsBootRefSam,     (*it)[i][j].transitions.fstRefSam()      );
	addReplicate( &transversionsBootRefSam,   (*it)[i][j].transversions.fstRefSam()    );
	addReplicate( &noDamageBootRefSam,        (*it)[i][j].noDamage.fstRefSam()         );

	addReplicate( &allBootSamRef,             (*it)[i][j].all.fstSamRef()           );
	addReplicate( &nocpgBootSamRef,           (*it)[i][j].noCpg.fstSamRef()         );
	addReplicate( &onlycpgBootSamRef,         (*it)[i][j].onlyCpg.fstSamRef()       );
	addReplicate( &transitionsBootSamRef,     (*it)[i][j].transitions.fstSamRef()   );
	addReplicate( &transversionsBootSamRef,   (*it)[i][j].transversions.fstSamRef() );
	addReplicate( &noDamageBootSamRef,        (*it)[i][j].noDamage.fstSamRef()      );

					    
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
	<<allDEVRefSam.first                <<"\t"<<allDEVRefSam.second            <<"\t"<<(all.fstRefSam()            / allDEVRefSam.second) <<"\t"
	<<allDEVSamRef.first                <<"\t"<<allDEVSamRef.second            <<"\t"<<(all.fstSamRef()            / allDEVSamRef.second) <<"\t"

	<<onlyCpg<<"\t"  
	<<onlyCpgDEVRefSam.first            <<"\t"<<onlyCpgDEVRefSam.second        <<"\t"<<(onlyCpg.fstRefSam()        / onlyCpgDEVRefSam.second) <<"\t"
	<<onlyCpgDEVSamRef.first            <<"\t"<<onlyCpgDEVSamRef.second        <<"\t"<<(onlyCpg.fstSamRef()        / onlyCpgDEVSamRef.second) <<"\t"

	<<noCpg<<"\t"  
	<<noCpgDEVRefSam.first              <<"\t"<<noCpgDEVRefSam.second          <<"\t"<<(noCpg.fstRefSam()          / noCpgDEVRefSam.second) <<"\t"
	<<noCpgDEVSamRef.first              <<"\t"<<noCpgDEVSamRef.second          <<"\t"<<(noCpg.fstSamRef()          / noCpgDEVSamRef.second) <<"\t"

	<<transitions<<"\t"  
	<<transitionsDEVRefSam.first        <<"\t"<<transitionsDEVRefSam.second    <<"\t"<<(transitions.fstRefSam()    / transitionsDEVRefSam.second) <<"\t"
	<<transitionsDEVSamRef.first        <<"\t"<<transitionsDEVSamRef.second    <<"\t"<<(transitions.fstSamRef()    / transitionsDEVSamRef.second) <<"\t"

	<<transversions<<"\t"  
	<<transversionsDEVRefSam.first      <<"\t"<<transversionsDEVRefSam.second  <<"\t"<<(transversions.fstRefSam()  / transversionsDEVRefSam.second) <<"\t"
	<<transversionsDEVSamRef.first      <<"\t"<<transversionsDEVSamRef.second  <<"\t"<<(transversions.fstSamRef()  / transversionsDEVSamRef.second) <<"\t"

	<<noDamage<<"\t"  
	<<noDamageDEVRefSam.first           <<"\t"<<noDamageDEVRefSam.second       <<"\t"<<(noDamage.fstRefSam()       / noDamageDEVRefSam.second) <<"\t"
	<<noDamageDEVSamRef.first           <<"\t"<<noDamageDEVSamRef.second       <<"\t"<<(noDamage.fstSamRef()       / noDamageDEVSamRef.second) ;
    

    return toreturn.str();
    
}
