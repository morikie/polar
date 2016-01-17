#!/bin/bash

grep "^.._" gbStatus.txt | cut -f 1,2 > ucsc_txRefSeq.txt 

