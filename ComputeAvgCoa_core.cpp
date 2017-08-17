/*
 * ComputeAvgCoa_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "ComputeAvgCoa_core.h"

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



inline bool isTransition(char allel_sample,char allel_reference,char allel_chimpHumanAncestor){
    return (  
	    (allel_sample == 'C' && (allel_reference == 'T' || (allel_chimpHumanAncestor == 'T') ) ) || 
	    (allel_sample == 'T' && (allel_reference == 'C' || (allel_chimpHumanAncestor == 'C') ) ) || 
	    (allel_sample == 'A' && (allel_reference == 'G' || (allel_chimpHumanAncestor == 'G') ) ) || 
	    (allel_sample == 'G' && (allel_reference == 'A' || (allel_chimpHumanAncestor == 'A') ) ) 
	      );
}


inline bool     isDamage(char allel_sample,char allel_reference,char allel_chimpHumanAncestor){
    return (  
	    (allel_sample == 'T' && (allel_reference == 'C' || (allel_chimpHumanAncestor == 'C') ) ) || 
	    (allel_sample == 'A' && (allel_reference == 'G' || (allel_chimpHumanAncestor == 'G') ) ) 		  
	      );
}

void computeDiv(char allel_chimpHumanAncestor,char allel_reference,char allel_sample,bool isCpG,AvgCoaResult * divr){
    if( (allel_chimpHumanAncestor == allel_reference ) && //no mutation
	(allel_chimpHumanAncestor == allel_sample ) ){
		
	//no need to call isTransition() isDamage() since no mutation
	divr->all.counterSame++;
	if(isCpG)
	    divr->onlyCpg.counterSame++;
	else
	    divr->noCpg.counterSame++;
    }


    if( (allel_chimpHumanAncestor != allel_reference ) || //mutation occured 
	(allel_chimpHumanAncestor != allel_sample ) ){



	if(allel_reference == allel_sample){ //common		    
#ifdef DEBUG3
	    cout<<"common"<<endl;
#endif


	    if(isTransition(allel_sample,allel_reference,allel_chimpHumanAncestor)){
		divr->transitions.counterCommon++;
	    }else{
		divr->transversions.counterCommon++;
	    }

	    if(!isDamage(    allel_sample,allel_reference,allel_chimpHumanAncestor)){
		divr->noDamage.counterCommon++;
	    }else{
	    }



	    divr->all.counterCommon++;
	    if(isCpG)
		divr->onlyCpg.counterCommon++;
	    else
		divr->noCpg.counterCommon++;

	}


	if( (allel_chimpHumanAncestor == allel_reference ) && //mutation occured in sample
	    (allel_chimpHumanAncestor != allel_sample ) ){

			    
#ifdef DEBUG3
	    cout<<"sample"<<endl;
#endif

	    if(isTransition(allel_sample,allel_reference,allel_chimpHumanAncestor)){
		divr->transitions.counterSample++;
	    }else{
		divr->transversions.counterSample++;

	    }

	    if(!isDamage(    allel_sample,allel_reference,allel_chimpHumanAncestor)){
		divr->noDamage.counterSample++;
	    }


	    divr->all.counterSample++;
	    if(isCpG)
		divr->onlyCpg.counterSample++;
	    else
		divr->noCpg.counterSample++;

	}

	if( (allel_chimpHumanAncestor != allel_reference ) && //mutation occured in reference
	    (allel_chimpHumanAncestor == allel_sample ) ){


			    
#ifdef DEBUG3
	    cout<<"reference"<<endl;
#endif



	    if(isTransition(allel_sample,allel_reference,allel_chimpHumanAncestor)){
		divr->transitions.counterReference++;
	    }else{
		divr->transversions.counterReference++;
	    }

	    if(!isDamage(    allel_sample,allel_reference,allel_chimpHumanAncestor)){
		divr->noDamage.counterReference++;
	    }



	    divr->all.counterReference++;
	    if(isCpG)
		divr->onlyCpg.counterReference++;
	    else
		divr->noCpg.counterReference++;

	}
				 
    }
			     


}
