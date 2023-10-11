"""
author: adam h freedman
afreedman405 at gmail.com
data: Fri Aug 26 10:55:18 EDT 2016

This script takes as an input Rcorrector error corrected Illumina paired-reads
in fastq format and:

1. Removes any reads that Rcorrector indentifes as containing an error,
but can't be corrected, typically low complexity sequences. For these,
the header contains 'unfixable'.

2. Strips the ' cor' from headers of reads that Rcorrector fixed, to avoid
issues created by certain header formats for downstream tools.

3. Write a log with counts of (a) read pairs that were removed because one end
was unfixable, (b) corrected left and right reads, (c) total number of
read pairs containing at least one corrected read.

Currently, this script only handles paired-end data, and handle either unzipped
or gzipped files on the fly, so long as the gzipped files end with 'gz'.

"""
import argparse
import sys        
import gzip

try:
    from itertools import izip_longest
except ImportError: 
    from itertools import zip_longest as izip_longest

import argparse
import os.path

def get_input_streams(r1file,r2file):
    if r1file[-2:]=='gz':
        r1handle=gzip.open(r1file,'rb')
        r2handle=gzip.open(r2file,'rb')
    else:
        r1handle=open(r1file,'r')
        r2handle=open(r2file,'r')
    
    return r1handle,r2handle
        
        
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)  
    

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for filtering and logging rCorrector fastq outputs")
    parser.add_argument('-1','--left_reads',dest='leftreads',type=str,help='R1 fastq file')
    parser.add_argument('-2','--right_reads',dest='rightreads',type=str,help='R2 fastq file')
    parser.add_argument('-s','--sample_id',dest='id',type=str,help='sample name to write to log file')
    opts = parser.parse_args()

    output_dir = os.path.dirname(opts.leftreads)
    r1out = open(os.path.join(output_dir, "unfixrm_{}".format(
        os.path.basename(opts.leftreads).replace('.gz', ''))), 'w')
    r2out = open(os.path.join(output_dir, "unfixrm_{}".format(
        os.path.basename(opts.rightreads).replace('.gz', ''))), 'w')
    r1_cor_count=0
    r2_cor_count=0
    pair_cor_count=0
    unfix_r1_count=0
    unfix_r2_count=0
    unfix_both_count=0   

    r1_stream,r2_stream=get_input_streams(opts.leftreads,opts.rightreads)

    with r1_stream as f1, r2_stream as f2:
        R1=grouper(f1,4)
        R2=grouper(f2,4)
        counter=0
        for entry in R1:
            counter+=1
            if counter%100000==0:
                print("%s reads processed" % counter)
        
            head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
            head2,seq2,placeholder2,qual2=[j.strip() for j in next(R2)]
            
            if 'unfixable' in head1 and 'unfixable' not in head2:
                unfix_r1_count+=1
            elif 'unfixable' in head2 and 'unfixable' not in head1:
                unfix_r2_count+=1
            elif 'unfixable' in head1 and 'unfixable' in head2:
                unfix_both_count+=1
            else:
                if 'cor' in head1:
                    r1_cor_count+=1
                if 'cor' in head2:
                    r2_cor_count+=1
                if 'cor' in head1 or 'cor' in head2:
                    pair_cor_count+=1
                
                r1out.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))
                r2out.write('%s\n' % '\n'.join([head2,seq2,placeholder2,qual2]))
    
    total_unfixable = unfix_r1_count+unfix_r2_count+unfix_both_count
    total_retained = counter - total_unfixable

    unfix_log = open(os.path.join(
        output_dir, "rmunfixable_{}.log".format(opts.id)), 'w')
    unfix_log.write('total PE reads:%s\nremoved PE reads:%s\nretained PE reads:%s\nR1 corrected:%s\nR2 corrected:%s\npairs corrected:%s\nR1 unfixable:%s\nR2 unfixable:%s\nboth reads unfixable:%s\n' % (counter,total_unfixable,total_retained,r1_cor_count,r2_cor_count,pair_cor_count,unfix_r1_count,unfix_r2_count,unfix_both_count))
            
    r1out.close()
    r2out.close() 
    unfix_log.close()
