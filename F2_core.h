/*
 * F2_core
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef F2_core_h
#define F2_core_h
#include <ctype.h>
#include "F2Result.h"
#include "libgab.h"

using namespace std;

//void computeF2(char allel_chimpHumanAncestor,char allel_reference,char allel_sample, bool isCpg, AvgCoaResult * divr);
//bool computeF2(const double freq_condition,const double freq_ind1,const double freq_ind2,const bool isCpG,F2Result * divr);
bool computeF2pair(const double freq_ind1,
		   const double freq_ind2,
		   const bool isCpG,
		   const bool isSitePotentialTransition, 
		   const bool isSitePotentialDamage, 
		   F2Result * f2res);

bool addCountersInd(const int refCount1,
		    const int altCount1,
		    const double freq_ind1,
		    const int refCount2,
		    const int altCount2,
		    const double freq_ind2,
		    const bool isCpG,
		    const bool isSitePotentialTransition, 
		    const bool isSitePotentialDamage,
		    F2Result * f2res);


#endif
