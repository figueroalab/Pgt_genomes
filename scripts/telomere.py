#!/usr/bin/env python

import sys
from Bio import SeqIO

def search_telomere( seq ):
    
    if seq.find("CCCTAACCCTAACCCTAA") >= 0:
        return 0-seq.count("CCCTAA")
    elif seq.find("TTAGGGTTAGGGTTAGGG" )>= 0:
        return seq.count("TTAGGG")
    else: return 0        

def main():
    if len(sys.argv) < 2:
        print("usage: %s in_fasta" % sys.argv[0], file=sys.stderr)                   
        #print("usage", sys.argv[0], file=sys.stderr)
        sys.exit()
   
    inhandle = open( sys.argv[1], "r" )
    for record in SeqIO.parse( inhandle, "fasta"):
        cnt = search_telomere( record.seq.upper() )
        if abs(cnt) >= 10:
            if cnt>0: print(record.seq[-200:]) 
            else: print(record.seq[:200])
            print("%s: %d" % (record.id, abs(cnt))) 

main()
