/*
 * F3_core
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "F3_core.h"
#include "ComputeAvgCoa_core.h"


bool computeF3(const double freq_condition,const double freq_ind1,const double freq_ind2,
	       const bool isCpG,//true if CpG
	       const bool isSitePotentialTransition, //true if transition
	       const bool isSitePotentialDamage, //true if damage
	       F3Result * f3res){

    
    double f3=(freq_condition-freq_ind1)*(freq_condition-freq_ind2);

    f3res->all.f3Sum                 += f3;
    f3res->all.counterSites ++;

    if(isCpG){
	f3res->onlyCpg.f3Sum         += f3;
	f3res->onlyCpg.counterSites  ++;

    }else{
	f3res->noCpg.f3Sum           +=f3;
	f3res->noCpg.counterSites++;
    }

    if(  isSitePotentialTransition ){
	f3res->transitions.f3Sum     +=f3;
	f3res->transitions.counterSites++;

    }else{
	f3res->transversions.f3Sum   +=f3;
	f3res->transversions.counterSites++;
    }
    if( !isSitePotentialDamage ){
	f3res->noDamage.f3Sum        +=f3;
	f3res->noDamage.counterSites++;
    }//TODO check!


    return false;
        
}
