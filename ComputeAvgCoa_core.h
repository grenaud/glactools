/*
 * ComputeAvgCoa_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ComputeAvgCoa_core_h
#define ComputeAvgCoa_core_h

#include "AvgCoaResult.h"

using namespace std;



inline bool isTransition(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
inline bool     isDamage(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
void computeDiv(char allel_chimpHumanAncestor,char allel_reference,char allel_sample, bool isCpg, AvgCoaResult * divr);

/* class ComputeAvgCoa_core{ */
/* private: */
/* public: */
/* ComputeAvgCoa_core(); */
/* ~ComputeAvgCoa_core(); */
/* }; */

#endif
