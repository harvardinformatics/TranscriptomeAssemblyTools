"""
author: adam h freedman
data: 03/2016
this script takes as an input Rcorrector error corrected reads
in fastq format and removes any reads that Rcorrect indentifes
as containing an error, but can't be corrected, typically
low complexity sequences. for these, the header contains 'unfixable'

upcated on 2017-05-10 to take either gzipped or unzipped fastq files.
gzipped files must end in 'gz' for the script to work properly!

"""

import sys        
import gzip
from itertools import izip,izip_longest

def get_input_streams(r1file,r2file):
    if sys.argv[1][-2:]=='gz':
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
    
r1out=open("unfixrm"+sys.argv[1].replace('.gz',''),'w')
r2out=open("unfixrm"+sys.argv[2].replace('.gz','') ,'w')
r1_cor_count=0
r2_cor_count=0
pair_cor_count=0
unfix_count=0   

r1_stream,r2_stream=get_input_streams(sys.argv[1],sys.argv[2])

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
            
            if 'unfixable' in head1 or 'unfixable' in head2:
                unfix_count+=1
            else:
                if 'cor' in head1:
                    r1_cor_count+=1
                if 'cor' in head2:
                    r2_cor_count+=1
                if 'cor' in head1 or 'cor' in head2:
                    pair_cor_count+=1
                head1=head1.replace(' cor','')
                head2=head2.replace(' cor','')
                r1out.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))
                r2out.write('%s\n' % '\n'.join([head2,seq2,placeholder2,qual2]))

unfix_log=open('rmunfixable.log','w')
unfix_log.write('total PE reads:%s\nremoved PE reads:%s\nretained PE reads:%s\nR1 corrected:%s\nR2 corrected:%s\npairs corrected:%s\n' % (counter,unfix_count,counter-unfix_count,r1_cor_count,r2_cor_count,pair_cor_count))
            
r1out.close()
r2out.close() 
