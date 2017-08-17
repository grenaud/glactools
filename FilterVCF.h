/*
 * FilterVCF
 * Date: Aug-22-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef FilterVCF_h
#define FilterVCF_h

#include "SimpleVCF.h"
#include "SetVCFFilters.h"

//extern int rejectIndel=0;
/* extern int rejectREFValidREF=0; */
/* extern int rejectREFValidALT=0; */
/* extern int rejectLOWCOV_REF=0; */
/* extern int rejectMap20=0; */
/* extern int rejectLOWMQ=0; */
/* extern int rejectLOWQUAL=0; */
/* extern int rejectCloseIndels=0; */
/* extern int rejectSysERR=0; */
/* extern int rejectRM=0; */
/* extern int rejectREF_unknownGeno=0; */

bool passedFilters(SimpleVCF * smvcf,const SetVCFFilters * filtersToUse);
//,int minCovcutoff,int maxCovcutoff,double minMapabilitycutoff,int minMQcutoff,int minGQcutoff);
string rejectFiltersTally();

#endif
