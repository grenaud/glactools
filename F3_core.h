/*
 * F3_core
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef F3_core_h
#define F3_core_h
#include <ctype.h>
#include "F3Result.h"
#include "utils.h"

using namespace std;

//void computeF3(char allel_chimpHumanAncestor,char allel_reference,char allel_sample, bool isCpg, AvgCoaResult * divr);
//bool computeF3(const double freq_condition,const double freq_ind1,const double freq_ind2,const bool isCpG,F3Result * divr);

bool computeF3triple(const double freq_condition,const double freq_ind1,const double freq_ind2,const bool isCpG,const bool isSitePotentialTransition, const bool isSitePotentialDamage, F3Result * f3res);

bool addCountersInd(const int refCount,
		    const int altCount,
		    const bool isCpG,
		    const bool isSitePotentialTransition, 
		    const bool isSitePotentialDamage,
		    F3Result * f3res);


#endif
