# TranscriptomeAssemblyTools
A collection of scripts for pre-processing fastq files prior to *de novo* transcriptome assembly, as well as slurm scripts for running what we view as a reasonable best practice workflow, the steps of which are:

1. run *fastqc* on raw fastq reads to identify potential issues with Illumina sequencing libraries such as retained adapter, over-represented (often low complexity) sequences, greater than expected decrease in base quality with read cycle
1. perform kmer-based read corrections with with [rCorrector](https://github.com/mourisl/Rcorrector), see [Song and Florea 2015, Gigascience](https://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0089-y)
1. remove read pairs where at least one read has been flagged by *rCorrector* as containing an erroneous kmer, and where it was not possible to computationally correct the errors
1. remove read pairs where at least read contains an over-represented sequence
1. perform light quality trimming with [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
1. optionally, to map the remaining reads to the [SILVA rRNA database](https://www.arb-silva.de/), then filtering read pairs where at least one read aligns to an rRNA sequence
1. assembly reads with [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

From our years of experience troubleshooting and evaluating *de novo* transcriptome assemblies, we have identified a number of statistical issues regarding the robustness of downstream analyses based upon them. A summary of these issues can be found at [Freedman et al. 2020, *Molecular Ecology Resources*](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13156). Given these issues and the rapidly decreasing cost of generating a genome assembly, we suggest that during the study design phase of your project, you consider the feasibility of assembling and annotating a genome before choosing to generate a *de novo* transcriptome assembly.

## Workflow steps
To the extent possible, to improve reproducibility (not to mention ease of implementation!), we run analyses from within conda environments. Below, we explain how to execute particular steps assuming that a separate conda enviornment is created for each step, with job scripts designed to be run on the SLURM job scheduler. These can easily be modified to work with other schedulers such as SGE and LSF.

### 1. Running fastqc
We can create a conda environment for *fastqc* as follows:
```bash
conda create -n fastqc -c bioconda fastqc
```
Then, we can run *fastqc*, one fastq file at a time, with [fastqc.sh](https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/job_scripts/fastqc.sh)
```bash
sbatch fastqc.sh <name of fastq file>
```
Note, if one supplies more threads to *fastqc*, with the *-t* switch, one can supply multiple fastq files at once. When choosing whether or not to run *fastqc* in multi-threaded mode, remember that the program only allocates 1 thread per file, i.e. there is no benefit to specifying more than one thread when only one file is being quality-checked. 


