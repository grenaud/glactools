#!/usr/bin/python

import sys;
import gzip;
from struct import *

#small example of how to read ACF/GLF files using python


DNAalphabet = ['N', 'A', 'C', 'G', 'T' ];

isACF=False;
isGLF=False;


fpz = gzip.open(sys.argv[1], 'rb');


#READING BINARY HEADER
magicnumber = fpz.read(4);
if(magicnumber != "BAM\1"):
    sys.stderr.write("The file does not begin with BAM\1\n");
    sys.exit(1);

formatcode = fpz.read(3);

isACF = (formatcode == "ACF");
isGLF = (formatcode == "GLF");

if(not(isACF) and not(isGLF)):
    sys.stderr.write("The file does not begin with ACF or GLF\n");
    sys.exit(1);

#bytesPerRecord = fpz.read(1);
bytesPerRecord = int(unpack('B', fpz.read(1) )[0]);
backslash1     = fpz.read(1);

sizeheader     = str(int(unpack('I', fpz.read(4) )[0]));


#UNPACKING TEXT HEADER
header         = fpz.read(int(sizeheader));
arrayChrName   = [];
for line in header.splitlines():
    print(line);
    values = line.split("\t")
    if(values[0] == "#SQ"):
        arrayChrName.append(values[1][3:]);



#READING SIZE POPULATIONS
sizepops       = int(unpack('I', fpz.read(4) )[0]);





#READING EACH "LINE"



try:
    while True:
        #COMMON FIELDS
        chriread      = fpz.read(2);
        if chriread == '':#we have reached the end of the file
            break;
        chri          = int(unpack('H', chriread )[0]);
        coord         = int(unpack('I', fpz.read(4) )[0]);
        refalt        = int(unpack('B', fpz.read(1) )[0]);
        ref           = refalt>>4;
        alt           = refalt & 15;

        strToPrint = "";
        #READ EACH INDIVIDUAL
        for p in range(0, sizepops+2):
            if( isACF ):
                acR  = int(unpack('H', fpz.read(bytesPerRecord) )[0]);
                acA  = int(unpack('H', fpz.read(bytesPerRecord) )[0]);
                cpg  = int(unpack('B', fpz.read(1) )[0]);
                strToPrint += str(acR)+","+str(acA)+":"+("1" if cpg else "0")+"\t";

            if( isGLF ):
                glRR = int(unpack('B', fpz.read(bytesPerRecord) )[0]);
                glRA = int(unpack('B', fpz.read(bytesPerRecord) )[0]);
                glAA = int(unpack('B', fpz.read(bytesPerRecord) )[0]);
                cpg  = int(unpack('B', fpz.read(1) )[0]);
                strToPrint += str(glRR)+","+str(glRA)+","+str(glAA)+":"+("1" if cpg else "0")+"\t";
        #PRINT DATA
        print arrayChrName[chri]+"\t"+str(coord)+"\t"+DNAalphabet[ref]+","+DNAalphabet[alt]+"\t"+strToPrint[:-1];
finally:
    fpz.close();
