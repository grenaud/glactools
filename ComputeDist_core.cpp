/*
 * ComputeDist_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "ComputeDist_core.h"

// #define DEBUG3


//          damage    (<-)
//        transition  (<->)
//    A    ---------   G  
//    |                |
//    |                |
//    |                |
//    |                |
//    C    ---------   T           
//       transition (<->)
//         damage   ( ->)



inline bool isTransition(char allel_1,char allel_2){
    return (  
	    (allel_1 == 'C' && allel_2 == 'T' ) || 
	    (allel_1 == 'T' && allel_2 == 'C' ) || 
	    (allel_1 == 'A' && allel_2 == 'G' ) || 
	    (allel_1 == 'G' && allel_2 == 'A' ) 
	      );
}


inline bool     isDamage(char allel_anc,char allel){
    //return isTransition(allel_1,allel_2);//we have no way of polarizing the direction of the mutation
    return (  
            (allel_anc == 'C' && allel == 'T' ) || 
            (allel_anc == 'G' && allel == 'A' )              
              );
}

void computeDist(const char allel_anc,const char allel_1,const char allel_2,bool isCpG,DistResult * divr){

    int allePairIndex=allelePair2Int(allel_1,allel_2);
    //cerr<<allel_1<<" "<<allel_2<<" "<<allePairIndex<<endl;
    divr->all.addAllelePair(allePairIndex);
    if(isCpG)
	divr->onlyCpg.addAllelePair(allePairIndex);
    else
	divr->noCpg.addAllelePair(allePairIndex);

    if(isTransition(allel_1,allel_2)){
	divr->transitions.addAllelePair(allePairIndex);
    }else{
	divr->transversions.addAllelePair(allePairIndex);
    }
    
    
    if( (allel_anc != 'N')              && 
	!isDamage( allel_anc,allel_1)   &&
	!isDamage( allel_anc,allel_2)  ){
	divr->noDamage.addAllelePair(allePairIndex);
    }else{
    }
			     
}
