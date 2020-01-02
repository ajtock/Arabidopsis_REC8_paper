#!/usr/bin/env python

# Author: Devon Ryan
# Description: a.k.a. foo.py as posted at https://www.biostars.org/p/95929/
#  Solves the problem of orphan reads in a pair that may remain after retaining
#  only those that align uniquely (i.e., after filtering using grep -v "XS:i:")
#  This script is invoked by samtools_PE_RNAseq_mapping_pipeline.sh

import csv
import sys

# ajt: added line for generating Unix (rather than Windows) text file
csv.register_dialect("unixEOL", delimiter="\t", lineterminator="\n")

f = csv.reader(sys.stdin, dialect="unixEOL")
of = csv.writer(sys.stdout, dialect="unixEOL")
last_read = None
for line in f :
    #take care of the header
    if(line[0][0] == "@") :
        of.writerow(line)
        continue

    if(last_read == None) :
        last_read = line
    else :
        if(last_read[0] == line[0]) :
            of.writerow(last_read)
            of.writerow(line)
            last_read = None
        else :
            last_read = line
