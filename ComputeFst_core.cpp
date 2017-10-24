/*
 * ComputeFst_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "ComputeFst_core.h"

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

void computeFst(char allel_chimpHumanAncestor,char allel_reference,char allel_sample,bool isCpG,FstResult * fstr){
    if( (allel_chimpHumanAncestor == allel_reference ) && //no mutation
	(allel_chimpHumanAncestor == allel_sample ) ){
		
	//no need to call isTransition() isDamage() since no mutation
	fstr->all.counterSame++;
	if(isCpG)
	    fstr->onlyCpg.counterSame++;
	else
	    fstr->noCpg.counterSame++;
    }


    if( (allel_chimpHumanAncestor != allel_reference ) || //mutation occured 
	(allel_chimpHumanAncestor != allel_sample ) ){



	if(allel_reference == allel_sample){ //common		    
#ifdef DEBUG3
	    cout<<"common"<<endl;
#endif


	    if(isTransition(allel_sample,allel_reference,allel_chimpHumanAncestor)){
		fstr->transitions.counterCommon++;
	    }else{
		fstr->transversions.counterCommon++;
	    }

	    if(!isDamage(    allel_sample,allel_reference,allel_chimpHumanAncestor)){
		fstr->noDamage.counterCommon++;
	    }else{
	    }



	    fstr->all.counterCommon++;
	    if(isCpG)
		fstr->onlyCpg.counterCommon++;
	    else
		fstr->noCpg.counterCommon++;

	}


	if( (allel_chimpHumanAncestor == allel_reference ) && //mutation occured in sample
	    (allel_chimpHumanAncestor != allel_sample ) ){

			    
#ifdef DEBUG3
	    cout<<"sample"<<endl;
#endif

	    if(isTransition(allel_sample,allel_reference,allel_chimpHumanAncestor)){
		fstr->transitions.counterSample++;
	    }else{
		fstr->transversions.counterSample++;

	    }

	    if(!isDamage(    allel_sample,allel_reference,allel_chimpHumanAncestor)){
		fstr->noDamage.counterSample++;
	    }


	    fstr->all.counterSample++;
	    if(isCpG)
		fstr->onlyCpg.counterSample++;
	    else
		fstr->noCpg.counterSample++;

	}

	if( (allel_chimpHumanAncestor != allel_reference ) && //mutation occured in reference
	    (allel_chimpHumanAncestor == allel_sample ) ){


			    
#ifdef DEBUG3
	    cout<<"reference"<<endl;
#endif



	    if(isTransition(allel_sample,allel_reference,allel_chimpHumanAncestor)){
		fstr->transitions.counterReference++;
	    }else{
		fstr->transversions.counterReference++;
	    }

	    if(!isDamage(    allel_sample,allel_reference,allel_chimpHumanAncestor)){
		fstr->noDamage.counterReference++;
	    }



	    fstr->all.counterReference++;
	    if(isCpG)
		fstr->onlyCpg.counterReference++;
	    else
		fstr->noCpg.counterReference++;

	}
				 
    }
			     


}
