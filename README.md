# TranscriptomeAssemblyTools
A collection of scripts for processing fastq files in ways to improve de novo transcriptome assemblies, and for evaluating those assemblies.

## FilterUncorrectablePEfastq.py
Takes a paired-end Illumina fastq file generated from an RNA-seq library, and that has been error corrected with rCorrector,and removes reads with errors but that are unfisable, and strips 'cor' flags from reads that were corrected.
