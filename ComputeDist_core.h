/*
 * ComputeDist_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ComputeDist_core_h
#define ComputeDist_core_h

#include "DistResult.h"

using namespace std;



inline bool isTransition(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
inline bool     isDamage(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
void computeDiv(char allel_chimpHumanAncestor,char allel_reference,char allel_sample, bool isCpg, DistResult * divr);

/* class ComputeDist_core{ */
/* private: */
/* public: */
/* ComputeDist_core(); */
/* ~ComputeDist_core(); */
/* }; */

#endif
