import sys        
import gzip
from itertools import izip,izip_longest
import argparse
from os.path import basename

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
    opts = parser.parse_args()

    r1out=open('sraheaderfixed_%s' % basename(opts.leftreads).replace('.gz',''),'w')
    r2out=open('sraheaderfixed_%s' % basename(opts.rightreads).replace('.gz',''),'w') 
    r1_stream,r2_stream=get_input_streams(opts.leftreads,opts.rightreads)
    with r1_stream as f1, r2_stream as f2:
        R1=grouper(f1,4)
        R2=grouper(f2,4)
        counter=0
        for entry in R1:
            counter+=1
            if counter%100000==0:
                print "%s reads processed" % counter
        
            head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
            head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()]
            head1=head1.split()[0]+'/1'
            head2=head2.split()[0]+'/2'
            placeholder1 = '+'
            placeholder2 = '+'
            r1out.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))
            r2out.write('%s\n' % '\n'.join([head2,seq2,placeholder2,qual2]))
            
    r1out.close()
    r2out.close() 
