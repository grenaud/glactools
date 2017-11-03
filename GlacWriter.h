/*
 * GlacWriter
 * Date: Jul-27-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef GlacWriter_h
#define GlacWriter_h

//#define DEBUGWRITEREFALT
//#define DEBUGWRITEAR
#include <cinttypes>
#include <string.h>

#include "htslib/bgzf.h"

#include "AlleleRecords.h"

using namespace std;

class GlacWriter{
private:
    uint32_t sizePops;
    bool glFormat;
    int bytesForRecord;
    bool uncompressed;
    BGZF * fpBGZF    ;
public:
    GlacWriter(uint32_t sizePops_,bool glFormat_,int bytesForRecord,int compressionThreads=1,bool uncompressed_=false);
    GlacWriter(const GlacWriter & other);
    ~GlacWriter();
    GlacWriter & operator= (const GlacWriter & other);
    bool writeHeader(string const & header) const;
    bool writeAlleleRecord(const AlleleRecords * toWrite) const;
};
#endif
