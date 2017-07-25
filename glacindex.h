/*
 * glacindex
 * Date: Jul-25-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef glacindex_h
#define glacindex_h

using namespace std;

#include <string>

class glacindex{
private:

public:
    glacindex();
    glacindex(const glacindex & other);
    ~glacindex();
    glacindex & operator= (const glacindex & other);
    
    string usage() const;
    int run(int argc, char *argv[]);

};
#endif
