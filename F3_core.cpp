/*
 * F3_core
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "F3_core.h"
#include "ComputeAvgCoa_core.h"


bool computeF3(const double freq_condition,const double freq_ind1,const double freq_ind2,const bool isCpG,F3Result * dSres){
    //must be different
    // if(allel_chimpHumanAncestor == allel_condition){
    // 	return false;
    // }

    // //require the reference allele to be either derived or ancestral
    // if( ( allel_ind1 == allel_chimpHumanAncestor) ||
    // 	( allel_ind1 == allel_condition)   ){
    // 	//fine
    // }else{
    // 	return false;
    // }

    // //require the sample allele to be either derived or ancestral
    // if( ( allel_ind2 == allel_chimpHumanAncestor) ||
    // 	( allel_ind2 == allel_condition)   ){
    // 	//fine
    // }else{
    // 	return false;
    // }

    //TODO

    double f3=(freq_condition-freq_ind1)*(freq_condition-freq_ind2);
    return false;
    
    //AA
    // if( ( allel_ind1    == allel_chimpHumanAncestor) &&
    // 	( allel_ind2    == allel_chimpHumanAncestor) ){


    // 	dSres->all.counterAncAnc++;
    // 	if(isCpG)
    // 	    dSres->onlyCpg.counterAncAnc++;
    // 	else
    // 	    dSres->noCpg.counterAncAnc++;

    // 	if(isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
    // 	    dSres->transitions.counterAncAnc++;
    // 	else
    // 	    dSres->transversions.counterAncAnc++;
	
    // 	if( (toupper(allel_chimpHumanAncestor)  == 'C' && 
    // 	     toupper(allel_condition)           == 'T'  ) 
    // 	    ||
    // 	    (toupper(allel_chimpHumanAncestor)  == 'G' && 
    // 	     toupper(allel_condition)           == 'A'  ) ){
    // 	    dSres->noDamage.counterAncAnc++;
    // 	}

    // 	//cout<<"AA"<<endl;
	
    // 	return false;

    // }else{
    // 	//AD
    // 	if( ( allel_ind1   ==  allel_chimpHumanAncestor) &&
    // 	    ( allel_ind2   ==  allel_condition) ){


    // 	    // cout<<"AD"<<endl;
    // 	    dSres->all.counterAncDer++;
    // 	    if(isCpG)
    // 		dSres->onlyCpg.counterAncDer++;
    // 	    else
    // 		dSres->noCpg.counterAncDer++;

    // 	    if(isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
    // 		dSres->transitions.counterAncDer++;
    // 	    else
    // 		dSres->transversions.counterAncDer++;

    // 	    if( (toupper(allel_chimpHumanAncestor)  == 'C' && 
    // 		 toupper(allel_condition)           == 'T'  ) 
    // 		||
    // 		(toupper(allel_chimpHumanAncestor)  == 'G' && 
    // 		 toupper(allel_condition)           == 'A'  ) ){
    // 		dSres->noDamage.counterAncDer++;
    // 	    }

    // 	    //cout<<"AD"<<endl;
	    
    // 	    return true;

    // 	}else{
    // 	    //DA
    // 	    if( ( allel_ind1  == allel_condition) &&
    // 		( allel_ind2  == allel_chimpHumanAncestor) ){

    // 		// cout<<"DA"<<endl;
    // 		dSres->all.counterDerAnc++;
    // 		if(isCpG)
    // 		    dSres->onlyCpg.counterDerAnc++;
    // 		else
    // 		    dSres->noCpg.counterDerAnc++;

    // 		if(isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
    // 		    dSres->transitions.counterDerAnc++;
    // 		else
    // 		    dSres->transversions.counterDerAnc++;

    // 		if( (toupper(allel_chimpHumanAncestor)  == 'C' && 
    // 		     toupper(allel_condition)           == 'T'  ) 
    // 		    ||
    // 		    (toupper(allel_chimpHumanAncestor)  == 'G' && 
    // 		     toupper(allel_condition)           == 'A'  ) ){
    // 		    dSres->noDamage.counterDerAnc++;
    // 		}

    // 		//cout<<"DA"<<endl;

    // 		return true;

    // 	    }else{

    // 		//DD
    // 		if( ( allel_ind1 == allel_condition) &&
    // 		    ( allel_ind2    == allel_condition) ){



    // 		    dSres->all.counterDerDer++;
    // 		    if(isCpG)
    // 			dSres->onlyCpg.counterDerDer++;
    // 		    else
    // 			dSres->noCpg.counterDerDer++;

    // 		    if(isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
    // 			dSres->transitions.counterDerDer++;
    // 		    else
    // 			dSres->transversions.counterDerDer++;


    // 		    if( (toupper(allel_chimpHumanAncestor)  == 'C' && 
    // 			 toupper(allel_condition)           == 'T'  ) 
    // 			||
    // 			(toupper(allel_chimpHumanAncestor)  == 'G' && 
    // 			 toupper(allel_condition)           == 'A'  ) ){
    // 			dSres->noDamage.counterDerDer++;
    // 		    }

    // 		    //cout<<"DD"<<endl;
		    
    // 		    return false;

    // 		}else{
    // 		    cerr<<"F3_core: Invalid state"<<endl;	     
    // 		    exit(1);
    // 		}


    // 	    }

    // 	}


    // }
    
}
