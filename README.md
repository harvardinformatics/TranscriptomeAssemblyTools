# TranscriptomeAssemblyTools
A collection of scripts for processing fastq files in ways to improve de novo transcriptome assemblies, and for evaluating those assemblies.

## FilterUncorrectablePEfastq.py
Takes a paired-end Illumina fastq file generated from an RNA-seq library, and that has been error corrected with rCorrector,and removes reads with errors but that are unfisable, and strips 'cor' flags from reads that were corrected.

## RemoveFastqcOverrepSequenceReads.py
Parses the fastqc output files to retrieve over-represented sequences, and uses these to remove read pairs where either read has a sequence match to an over-represented sequence.
