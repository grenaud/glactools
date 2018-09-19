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



inline bool isTransition(const char allel_1,const char allel_2);
inline bool     isDamage(const char allel_1,const char allel_2);
void computeDist(const char allel_anc,const char allel_1,const char allel_2,bool isCpg, DistResult * divr);
void addIdenticalSite(const char allel_anc,const char allel_1,bool isCpG,DistResult * divr,int allePairIndex);
/* class ComputeDist_core{ */
/* private: */
/* public: */
/* ComputeDist_core(); */
/* ~ComputeDist_core(); */
/* }; */

#endif
