# TranscriptomeAssemblyTools
A collection of scripts for pre-processing fastq files prior to *de novo* transcriptome assembly, as well as slurm scripts for running what we view as a reasonable best practice workflow, the steps of which are:

1. run *fastqc* on raw fastq reads to identify potential issues with Illumina sequencing libraries such as retained adapter, over-represented (often low complexity) sequences, greater than expected decrease in base quality with read cycle
1. perform kmer-based read corrections with with [rCorrector](https://github.com/mourisl/Rcorrector), see [Song and Florea 2015, Gigascience](https://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0089-y)
1. remove read pairs where at least one read has been flagged by *rcorrector* as containing an erroneous kmer, and where it was not possible to computationally correct the errors
1. remove read pairs where at least read contains an over-represented sequence
1. perform light quality trimming with *trimgalore*
1. optionally, to map the remaining reads to the SILVA rRNA database, then filtering read pairs where at least one read aligns to an rRNA sequence
1. assembly reads with Trinity

