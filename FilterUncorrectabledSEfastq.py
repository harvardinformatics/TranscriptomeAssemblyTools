"""
author: Benedek Dankó
2020-11-21

This script takes as an input Rcorrector error corrected Illumina single-reads
in fastq format and:

1. Removes any reads that Rcorrector indentifes as containing an error,
but can't be corrected, typically low complexity sequences. For these,
the header contains 'unfixable'.

2. Strips the ' cor' from headers of reads that Rcorrector fixed, to avoid
issues created by certain header formats for downstream tools.

3. Writes a log.

Currently, this script handles single-end data, and handle either unzipped
or gzipped files on the fly, so long as the gzipped files end with 'gz'.

"""

import sys        
import gzip
from itertools import izip,izip_longest
import argparse
from os.path import basename

def get_input_streams(file):
    if file[-2:]=='gz':
        handle=gzip.open(file,'rb')
    else:
        handle=open(file,'r')
    return handle
        
        
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)  
    

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for filtering and logging rCorrector fastq outputs")
    parser.add_argument('-f','--reads', dest='reads',type=str,help='SE fastq file')
    parser.add_argument('-s','--sample_id', dest='id',type=str,help='sample name to write to log file')
    opts = parser.parse_args()

    out=open('%s' % basename(opts.reads).replace('.gz',''),'w')

    cor_count=0
    unfix_count=0 

    stream=get_input_streams(opts.reads)

    with stream as f:
        R=grouper(f,4)
        counter=0
        for entry in R:
            counter+=1
            if counter%100000==0:
                print ("%s reads processed" % counter)
        
            head,seq,placeholder,qual=[i.strip() for i in entry]
            
            if 'unfixable' in head:
                unfix_count+=1
            else:
                if 'cor' in head:
                    cor_count+=1
                
                head=head.split('l:')[0][:-1] 
                out.write('%s\n' % '\n'.join([head,seq,placeholder,qual]))
    
    total_unfixable = unfix_count
    total_retained = counter - total_unfixable

    unfix_log=open('rmunfixable_%s.log' % opts.id,'w')
    unfix_log.write('Total reads:%s\nRemoved reads:%s\nRetained reads:%s\nCorrected:%s\nUnfixable:%s\n' % (counter,total_unfixable,total_retained,cor_count,unfix_count))
            
    out.close()
    unfix_log.close()
