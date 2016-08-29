import sys
import gzip
from os.path import basename
import argparse
import re
from itertools import izip,izip_longest


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


def ParseKmerLog(fastqclog):    
    kmer_dict={}
    with open(fastqclog) as fp:
        for result in re.findall('Kmer Content(.*?)END_MODULE', fp.read(), re.S):
            for kmer in result.split('\n')[2:-1]:
                kmer_dict[kmer.split('\t')[0]]={}
                kmer_dict[kmer.split('\t')[0]]['count']=int(kmer.split('\t')[1])
                kmer_dict[kmer.split('\t')[0]]['P']=float(kmer.split('\t')[2])
                kmer_dict[kmer.split('\t')[0]]['obsexpmax']=float(kmer.split('\t')[3])
                kmer_dict[kmer.split('\t')[0]]['position']=float(kmer.split('\t')[4])  
                kmer_dict[kmer.split('\t')[0]]['removecount']=round(kmer_dict[kmer.split('\t')[0]]['count']-(kmer_dict[kmer.split('\t')[0]]['count']/float(kmer_dict[kmer.split('\t')[0]]['obsexpmax'i])))
    return kmer_dict


if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for kmer cleaning of fastq file to rm seq with overrep kmers")
    parser.add_argument('-1','--left_reads',dest='leftreads',type=str,help='R1 fastq file')
    parser.add_argument('-2','--right_reads',dest='rightreads',type=str,help='R2 fastq file')
    parser.add_argument('-fql','--fastqc_left',dest='l_fastqc',type=str,help='fastqc text file for R1')
    parser.add_argument('-fqr','--fastqc_right',dest='r_fastqc',type=str,help='fastqc text file for R2')
    parser.add_argument('-L','--seq_length',dest='seqlen',type=str,help='length of read sub-sequence to search for')
    parser.add_argument('-t','--count_threshold',dest='thresh',default=1,type=int,help='excess count threshold relative to expectation')
    opts = parser.parse_args() 
    
    left_kmer_seq_dict={}
    right_kmer_seq_dict={}

    leftseqs=ParseKmerLog(opts.l_fastqc)
                
    
    rightseqs=ParseKmerLog(opts.r_fastqc)
       
   
    r1_out=open('kmerclean_'+basename(opts.leftreads).replace('.gz',''),'w')
    r2_out=open('kmerclean_'+basename(opts.rightreads).replace('.gz',''),'w')
    
    r1_stream,r2_stream=get_input_streams(opts.leftreads,opts.rightreads)


    with r1_stream as f1, r2_stream as f2:
        R1=FastqIterate(f1)
        R2=FastqIterate(f2)
        for entry in R1:
            counter+=1
            if counter%100000==0:
                print "%s reads processed" % counter
        
            head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
            head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()] 

            
             
