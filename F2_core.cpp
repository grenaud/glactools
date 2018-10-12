/*
 * F2_core
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "F2_core.h"
#include "ComputeAvgCoa_core.h"




bool addCountersInd(const int  refCount1,
		    const int  altCount1,
		    const double freq_ind1,
		    const int  refCount2,
		    const int  altCount2,
		    const double freq_ind2,
		    const bool isCpG,
		    const bool isSitePotentialTransition, 
		    const bool isSitePotentialDamage,
		    F2Result * f2res){
    //cerr<<"addCountersInd "<<refCount<<" "<<altCount<<endl;
    f2res->all.addAlleleCounts(refCount1,altCount1,freq_ind1,refCount2,altCount2,freq_ind2);
    if(isCpG){
	f2res->onlyCpg.addAlleleCounts(refCount1,altCount1,freq_ind1,refCount2,altCount2,freq_ind2);
    }else{
	f2res->noCpg.addAlleleCounts(refCount1,altCount1,freq_ind1,refCount2,altCount2,freq_ind2);
    }

    if(  isSitePotentialTransition ){
	f2res->transitions.addAlleleCounts(refCount1,altCount1,freq_ind1,refCount2,altCount2,freq_ind2);
    }else{
	f2res->transversions.addAlleleCounts(refCount1,altCount1,freq_ind1,refCount2,altCount2,freq_ind2);
    }
    if( !isSitePotentialDamage ){
	f2res->noDamage.addAlleleCounts(refCount1,altCount1,freq_ind1,refCount2,altCount2,freq_ind2);
    }//TODO check!
        
}

bool computeF2pair(const double freq_ind1,
		   const double freq_ind2,
		   const bool isCpG,//true if CpG
		   const bool isSitePotentialTransition, //true if transition
		   const bool isSitePotentialDamage, //true if damage
		   F2Result * f2res){
    
    
    double f2=(freq_ind1-freq_ind2)*(freq_ind1-freq_ind2);//(f1-f2)^2
    //cerr<<f3<<"\t"<<freq_condition<<"\t"<<freq_ind1<<"\t"<<freq_condition<<"\t"<<freq_ind2<<"\t"<<f3res->all.f3Sum<<endl;

    f2res->all.f2Sum                 += f2;
    
    //cerr<<f3<<"\t"<<freq_condition<<"\t"<<freq_ind1<<"\t"<<freq_condition<<"\t"<<freq_ind2<<"\t"<<f3res->all.f3Sum<<endl;

    f2res->all.counterSites ++;

    if(isCpG){
	f2res->onlyCpg.f2Sum         += f2;
	f2res->onlyCpg.counterSites  ++;

    }else{
	f2res->noCpg.f2Sum           +=f2;
	f2res->noCpg.counterSites++;
    }

    if(  isSitePotentialTransition ){
	f2res->transitions.f2Sum     +=f2;
	f2res->transitions.counterSites++;

    }else{
	f2res->transversions.f2Sum   +=f2;
	f2res->transversions.counterSites++;
    }
    if( !isSitePotentialDamage ){
	f2res->noDamage.f2Sum        +=f2;
	f2res->noDamage.counterSites++;
    }//TODO check!


    return false;
        
}
