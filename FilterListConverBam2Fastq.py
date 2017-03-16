import argparse
from subprocess import Popen,PIPE
import glob
from sets import Set
from os import path


# set bits for reads in properly mapped pairs 
mate_to_bits={'R1':'0x42','R2':'0x82'}

def sortbam(unsorted_bam):
    """
    coordinate sorts bam file
    """
    print('coordinate sorting %s\n' % unsorted_bam)
    cmd='samtools sort %s> sorted_%s' % (unsorted_bam,unsorted_bam)
    samsort=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    samsort_stdout,samsort_stderr=samsort.communicate()
    if samsort.returncode==0:
        return True
    else:
        return False    


def extractbamheader(bamin):
    cmd='samtools view -H %s > header.sam' % bamin
    headgrab=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    headout,headerr=headgrab.communicate()
    if headgrab.returncode==0:
        headpass=True
        headerr=''
    else:
        headpass=False
    return headpass,headerr


def addheadertosam(samin):
    concatcmd='cat header.sam namesort_%s > wheader_%s' % (samin,samin)
    doconcat=Popen(concatcmd,shell=True,stderr=PIPE,stdout=PIPE)
    concatout,concaterr=doconcat.communicate()
    if doconcat.returncode==0:
        concatpass=True
    else:
        concatpass=False
        raise ValueError('%s\n' % concaterr)
    return concatpass

def sam2bam(samin):
    bamcmd='samtools view -Sbh %s > %s.bam' % (samin,samin[:-4])
    makebam=Popen(bamcmd,shell=True,stederr=PIPE,stdout=PIPE)
    bamout,bamerr=makebam.communicate()
    if makebam.returncode==0:
        makebampass=True
    else:
        makebampass=False
        raise ValueError('%s\n' % makebamerr)
    return makebampass
    

def namesortsam(samin):
    """
    takes filtered sam file, name sorts it,
    and writes it, with header, as a bam file
    """
    sortcmd='samtools sort -n -o namesort_%s %s' % (samin,samin)
    dosort=Popen(sortcmd,shell=True,stderr=PIPE,stdout=PIPE)
    sortout,sorterr=dosort.communicate()
    if dosort.returncode==0:
        sortpass=True
        sorterr=''
    else:
        sortpass=False
        #raise ValueError('%s\n' % sorterr)    
    return    sortpass,sorterr
 
def indexbam(sorted_bam):
    print('indexing %s\n' % sorted_bam) 
    cmd='samtools index %s' % sorted_bam
    samindex=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    samindex_stdout,samindex_stderr=samindex.communicate()
    if samindex.returncode==0:
        return True,''
    else:
        return True,samindex_stderr

def makefilterset(filterfile):
    fopen=open(filterfile,'r')
    filterset=Set()
    for line in fopen:
        contig=line.strip()
        filterset.add(contig)
    return filterset

def write_mate_sams(bamin,hexval,mate):
    """
    extract mapped reads from a coord-sorted
    bam file for a pass-filter contig in a
     mate-specific fashion
    """
    cmd='samtools view -h -f %s %s > %s_%s.sam' % (hexval,bamin,bamin[:-4],mate)
    getreads=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    readsout,readserr=getreads.communicate()
    if getreads.returncode==0:
        return True,''
    else:
        return False,readserr

def sam2fastq(samfile,contigskeep):
    counter=0
    fopen=open(samfile,'r')
    fqout=open('%s.fq' % samfile[:-4],'w')
    lastid=''
    for line in fopen:
        linelist=line.strip().split('\t')    
        if line[0]!='@' and linelist[2] in contigskeep:
            counter+=1
            if counter%100000==0:
                print 'processing read %s' % counter
            if lastid=='':
                lastid=linelist[0]
                fqout.write('%s\n' % '\n'.join([linelist[0],linelist[9],'+',linelist[10]]))
            elif lastid!=linelist[0]:
                lastid=linelist[0]
                fqout.write('%s\n' % '\n'.join([linelist[0],linelist[9],'+',linelist[10]]))
            else:
                pass
        else:
            pass

    fqout.close()
          

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="using contig keep list to extract fastq files")
    parser.add_argument('-s','--is_coord_sorted',action='store_true',help='flag for whether bam is name sorted')
    parser.add_argument('-b','--bam_in',dest='bamin',type=str,help='compressed read alignment infile')
    parser.add_argument('-k','--contig_keep_file',dest='keeps',type=str,help='list of contig mapping to keep')
    opts = parser.parse_args()  
    
    print opts
    
    ### make sure in bam is coord sorted to enable contig searhes   
    if opts.is_coord_sorted:
        sortedbam=opts.bamin
    else:
        bamsorter=sortbam(opts.bamin)
        sortedbam='sorted_%s' % opts.bamin
    
    ### verify or make bam file index is present
    if glob.glob('%s.bai' % sortedbam)==[]:
        indexflag,indexerr=indexbam(sortedbam) 
        if indexflag==False:
            raise ValueError('%s\n' % indexerr)

    ### make list of pass-filter contigs
    good_contigs=makefilterset(opts.keeps)

    ### write a contig-filtered bam file for each mate in a pair: 
    for mate in ['R1','R2']:
        print 'writing %s sam file' % mate
        mateflag,materr=write_mate_sams(sortedbam,mate_to_bits[mate],mate)
        if mateflag==True:
            pass
        else:
            raise ValueError('%s\n' % materr)

        print 'name-sorting %s sam file' % mate
        sortflag,sorterr=namesortsam('%s_%s.sam' % (sortedbam[:-4],mate))
        if sortflag==False:
            raise ValueError('%s\n' % sorterr)
        
        else:
            print 'writing unique %s reads to fastq' % mate
            sam2fastq('namesort_%s_%s.sam' % (sortedbam[:-4],mate),good_contigs)




