/*
 * Dstat_core
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Dstat_core_h
#define Dstat_core_h
#include <ctype.h>
#include "DstatResult.h"
#include "libgab.h"

using namespace std;

//void computeDstat(char allel_chimpHumanAncestor,char allel_reference,char allel_sample, bool isCpg, AvgCoaResult * divr);
bool computeDstat(const char allel_chimpHumanAncestor,const char allel_condition,const char allel_ind1,const char allel_ind2,const bool isCpG,DstatResult * divr);

#endif
