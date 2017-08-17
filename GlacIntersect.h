/*
 * GlacIntersect
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacIntersect_h
#define GlacIntersect_h

//#define DEBUG
#include <string>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"
#include "GlactoolsOperations.h"

using namespace std;


class GlacIntersect{
private:
    bool uncompressed=0;    
    int minPLdiffind=33;

public:
    GlacIntersect();
    GlacIntersect(const GlacIntersect & other);
    ~GlacIntersect();
    GlacIntersect & operator= (const GlacIntersect & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
