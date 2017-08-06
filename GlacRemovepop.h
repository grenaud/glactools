/*
 * GlacRemovepop
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacRemovepop_h
#define GlacRemovepop_h

#include <string>
#include <set>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GlacRemovepop{
private:
    bool uncompressed=0;    

public:
    GlacRemovepop();
    GlacRemovepop(const GlacRemovepop & other);
    ~GlacRemovepop();
    GlacRemovepop & operator= (const GlacRemovepop & other);
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
