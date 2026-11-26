/*
 * SafeJackknife
 *
 */

#ifndef SafeJackknife_h
#define SafeJackknife_h

#include <cmath>
#include <limits>
#include <vector>

#include "libgab.h"

using namespace std;

//Delete-one jackknife interval, or (nan,nan) when it is not defined: fewer than two blocks, or an
//estimate that cannot be computed (empty category) for all the data or for one of the replicates.
inline pair<double,double> jackknifeOrNan(double estimate,const vector<double> & replicates){
    bool defined = replicates.size() >= 2 && std::isfinite(estimate);
    for(unsigned int i=0;i<replicates.size();i++)
	defined = defined && std::isfinite(replicates[i]);

    if(!defined){
	const double notDefined = std::numeric_limits<double>::quiet_NaN();
	return pair<double,double>(notDefined,notDefined);
    }

    return computeJackknifeConfIntervals(estimate,replicates);
}

#endif
