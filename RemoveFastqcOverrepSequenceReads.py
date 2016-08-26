import sys
import gzip
from os.path import basename
import argparse
import re
from itertools import izip,izip_longest

def seqsmatch(overreplist,read):
    flag=False
    if overreplist!=[]:
        for seq in overreplist:
            if seq in read:
                flag=True
                break
    return flag
            
def get_input_streams(r1file,r2file):
    if  r1file[-2:]=='gz':
        r1handle=gzip.open(r1file,'rb')
        r2handle=gzip.open(r2file,'rb')
    else:
        r1handle=open(r1file,'r')
        r2handle=open(r2file,'r')
    
    return r1handle,r2handle

def FastqIterate(iterable,fillvalue=None):
    "Grab one 4-line fastq read at a time"
    args = [iter(iterable)] * 4
    return izip_longest(fillvalue=fillvalue, *args) 

def ParseFastqcLog(fastqclog):    
    with open(fastqclog) as fp:
        for result in re.findall('Overrepresented sequences(.*?)END_MODULE', fp.read(), re.S):
            seqs=([i.split('\t')[0] for i in result.split('\n')[2:-1]])
    return seqs     

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for removing reads with over-represented sequences")
    parser.add_argument('-1','--left_reads',dest='leftreads',type=str,help='R1 fastq file')
    parser.add_argument('-2','--right_reads',dest='rightreads',type=str,help='R2 fastq file')
    parser.add_argument('-fql','--fastqc_left',dest='l_fastqc',type=str,help='fastqc text file for R1')
    parser.add_argument('-fqr','--fastqc_right',dest='r_fastqc',type=str,help='fastqc text file for R2')
    opts = parser.parse_args()

    leftseqs=ParseFastqcLog(opts.l_fastqc)
    rightseqs=ParseFastqcLog(opts.r_fastqc)
   
    r1_out=open('rmoverrep_'+basename(opts.leftreads).replace('.gz',''),'w')
    r2_out=open('rmoverrep_'+basename(opts.rightreads).replace('.gz',''),'w')
    
    r1_stream,r2_stream=get_input_streams(opts.leftreads,opts.rightreads)
        
    counter=0
    failcounter=0

    with r1_stream as f1, r2_stream as f2:
        R1=FastqIterate(f1)
        R2=FastqIterate(f2)
        for entry in R1:
            counter+=1
            if counter%100000==0:
                print "%s reads processed" % counter
        
            head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
            head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()]
            
            flagleft,flagright=seqsmatch(leftseqs,seq1),seqsmatch(rightseqs,seq2)
            
            if True not in (flagleft,flagright):
                r1_out.write('%s\n' % '\n'.join([head1,seq1,'+',qual1]))
                r2_out.write('%s\n' % '\n'.join([head2,seq2,'+',qual2]))
            else:
                failcounter+=1


        print 'total # of reads evaluated = %s' % counter
        print 'number of reads retained = %s' % (counter-failcounter)
        print 'number of PE reads filtered = %s' % failcounter


r1_out.close()
r2_out.close()
