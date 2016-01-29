#!/bin/bash

#extracting first and second column from gbStatus where the lines start with "XX_" (X being any character, e.g. "NM_")
grep "^.._" gbStatus.txt | cut -f 1,2 > ucsc_txRefSeq.txt 

