/*
 * ComputeFst_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ComputeFst_core_h
#define ComputeFst_core_h

#include "FstResult.h"

using namespace std;



inline bool isTransition(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
inline bool     isDamage(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
void computeFst(char allel_chimpHumanAncestor,char allel_reference,char allel_sample,bool isCpG,FstResult * fstr);

/* class ComputeAvgCoa_core{ */
/* private: */
/* public: */
/* ComputeAvgCoa_core(); */
/* ~ComputeAvgCoa_core(); */
/* }; */

#endif
