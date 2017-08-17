/*
 * GLF2ACF
 * Date: Jul-30-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GLF2ACF_h
#define GLF2ACF_h

#include <string>

#include "utils.h"
#include "GlacWriter.h"
#include "GlacParser.h"
#include "AlleleRecords.h"

using namespace std;

class GLF2ACF{
private:
    bool uncompressed=0;    
    int minPLdiffind=33;

public:
    GLF2ACF();
    GLF2ACF(const GLF2ACF & other);
    ~GLF2ACF();
    GLF2ACF & operator= (const GLF2ACF & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
