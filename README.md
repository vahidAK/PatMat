# PatMat:  
![](docs/FlowChart.png)  

This workflow enables simultanous chromosome-scale haplotyping and parent-of-origin detection using a combination of nanopore sequencing and Strand-seq.
We will use nanopore-detected variants and their long-range phasing from Strand-seq to detect chromosome-scale haplotypes. We then use DNA methylation information at known imprinted regions to detect parent-of-origin.  

To run this workflow you will need the following third-party tools for processing nanopore data:  
**Guppy**: For basecalling nanopore reads.  
**[Minimap2](https://github.com/lh3/minimap2)**: To align nanopore reads to the reference genome.  
**[Nanopolish](https://github.com/jts/nanopolish)**: To call DNA methylation from nanopore data.  
**[NanoMethPhase](https://github.com/vahidAK/NanoMethPhase)**: To process methylation call results from nanopolish.  
**[Clair3](https://github.com/HKU-BAL/Clair3)**: To call variants from aligned nanopore reads.  

The workflow was developed using the above tools. However, you may use alternatives to each tool. Additional software tools are required to use Strand-seq (described in part 2).

Finally you need to use our tool in this repository "**PatMat.py**" to detect chromosome-scale parent-of-origin resolved haplotypes. You can clone the GitHub repository and use the PatMat.py. Before using PatMat.py you need to satisfy the following dependencies:  
[bgzip](http://www.htslib.org/doc/bgzip.html)  
[tabix](http://www.htslib.org/doc/tabix.html)  
python>=3.7.4 and its following dependencies:  
pytabix>=0.1  
pysam>=0.16.0  
tqdm>=4.54.1  

To run the Strand-seq part you need to satisfy dependencies described in the [scripts/Strand-seq/](https://github.com/vahidAK/PatMat/tree/main/scripts/Strand-seq) README.   

Table of Contents
=================

**[Full Tutorial](https://github.com/vahidAK/PatMat/blob/main/README.md#Full-Tutorial)**
* [Nanopore Data Analysis](https://github.com/vahidAK/PatMat/blob/main/README.md#1--Nanopore-Data-Analysis)
  * [Basecalling from nanopore data](https://github.com/vahidAK/PatMat/blob/main/README.md#1\-1--Basecalling-from-nanopore-data)
  * [Mapping nanopore basecalled reads](https://github.com/vahidAK/PatMat/blob/main/README.md#1\-2-Mapping-nanopore-basecalled-reads)
  * [Methylation Calling from nanopore data](https://github.com/vahidAK/PatMat/blob/main/README.md#1\-3-Methylation-Calling-from-nanopore-data)
  * [Variant Calling from nanopore data](https://github.com/vahidAK/PatMat/blob/main/README.md#1\-4-Variant-Calling-from-nanopore-data)
* [Strand-seq Data Analysis](https://github.com/vahidAK/PatMat/blob/main/README.md#2--Strand\-seq-Data-Analysis)
* [Parent-of-origin detection](https://github.com/vahidAK/PatMat/blob/main/README.md#3--Parent\-of\-origin-detection)
  
# Full Tutorial

## 1- Nanopore Data Analysis
### 1-1 Basecalling from nanopore data
We use Oxford Nanopore Technologies' basecaller, guppy, for translating raw signals to DNA sequence.
```
guppy_basecaller --input_path <Path to fast5 directory> \
  --save_path <path to output directory> \
  --config <appropriate configuration file (e.g. dna_r9.4.1_450bps_sup_prom.cfg)> \
  --device <GPU devices to be used (e,g. cuda:0 cuda:1)> \
  --trim_strategy dna
```
After basecalling you need to merge all the fastq files to a single file.

### 1-2 Mapping nanopore basecalled reads
Here we use [minimap2](https://github.com/lh3/minimap2) to align nanopore reads to reference genome. You may use other tools such as [Winnowmap](https://github.com/marbl/Winnowmap).
```
minimap2 -ax map-ont --MD -L -t <# of threads> \
  /path/to/reference.fa \
  /path/to/Nanopore_reads.fastq > /path/to/Nanopore_aligned_reads.sam 
```
After alignment was complete you need to sort and index alignment file.
```
samtools sort -@ <# of threads> /path/to/Nanopore_aligned_reads.sam  -o /path/to//Nanopore_aligned_reads.bam 
samtools index -@ <# of threads> /path/to/Nanopore_aligned_reads.bam
```

### 1-3 Methylation Calling from nanopore data  
Here we use [nanopolish](https://github.com/jts/nanopolish) for methylation calling however you may use [megalodon](https://github.com/nanoporetech/megalodon) or [DeepSignal](https://github.com/bioinfomaticsCSU/deepsignal). 

#### 1-3-1 indexing fastq file and fast5 files:

NOTE: Fastqs must be merged to a single file

```
nanopolish index -d /path/to/Nanopore_reads.fastq
```

#### 1-3-2 Methylation calling for CpG from each read:

```
nanopolish call-methylation \
  -t <number_of_threads> -q cpg \
  -r /path/to/Nanopore_reads.fastq \
  -b /path/to/Nanopore_aligned_reads.bam \
  -g /path/to/reference.fa > /path/to/MethylationCall.tsv
```

#### 1-3-3 Pre-processing methylation call file
We then need to pre-process methylation call file from nanopolish using [NanoMethPhase](https://github.com/vahidAK/NanoMethPhase) methyl_call_processor module.
```
nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > NanoMethPhase_MethylationCall.bed.gz && tabix -p bed NanoMethPhase_MethylationCall.bed.gz
```
### 1-4 Variant Calling from nanopore data

Here use [Clair3](https://github.com/HKU-BAL/Clair3) to call variants. However, you may call variants with other
tools such as [deepvariant](https://github.com/google/deepvariant).

```
run_clair3.sh --bam_fn=/path/to/Nanopore_aligned_reads.bam \
  --ref_fn=/path/to/reference.fa \
  --output=/path/to/output/directory \
  --threads=<# of threads> --platform=ont \
  --model_path=/path/to/model/ont_guppy5_r941_sup_g5014
```
After variant calling the results will be in merge_output.vcf.gz file in the output directory. You then need to extract high uality variants:  
```
gunzip -c /path/to/output/directory/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /path/to/output/Passed_Clair3_Variants.vcf
```  

## 2- Strand-seq Data Analysis
Typically, 30-100 good quality Strand-seq libraries with at least 20 million unique reads in total should be used for phasing. These must be aligned to the reference genome and poor-quality libraries must be identified and removed using [ASHLEYS QC](https://github.com/friendsofstrandseq/ashleys-qc) (Gros et al. 2021). Then, with a nanopore-derived VCF file of SNVs, the scripts in [scripts/Strand-seq/](https://github.com/vahidAK/PatMat/tree/main/scripts/Strand-seq) must be run in the directory containing the Strand-seq BAM files. This requires some installations (described in the [scripts/Strand-seq/](https://github.com/vahidAK/PatMat/tree/main/scripts/Strand-seq) README) as well as altering the header of the master.sh file. The command is then simply:
```
bash master.sh
```
This master script calls inversions using [InvertypeR](https://github.com/vincent-hanlon/InvertypeR) (most of the runtime), which help refine phasing, and then it phases the SNVs using the standard Strand-seq R packages [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) and [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR). The result is a phased VCF file of SNVs ("samplename.phased.inv_aware.vcf"), which can be used with PatMat.py as described below.

## 3- Parent-of-origin detection
Finally, parent-of-origin chromosome-scale haplotypes can be built using PatMat.py:  
```
PatMat.py phase -v /path/to/Passed_Clair3_Variants.vcf \
 -sv /path/to/StrandSeq_phased_variants.vcf \
 -mc /path/to/NanoMethPhase_MethylationCall.bed.gz \
 -b /path/to/Nanopore_aligned_reads.bam \
 -o <Output pass and prefix> \
 -t <# of threads>
```

Here is the full list of options:  
```
python patmat/PatMat.py phase -h

usage: nanomethphase phase --bam BAM --output OUTPUT --vcf VCF --strand_vcf
                           PHASED_VCF [--known_dmr KNOWN_DMR]
                           [--methylcallfile METHYLCALLFILE] [-h]
                           [--whatshap_vcf WHATSHAP_VCF]
                           [--whatshap_block WHATSHAP_BLOCK]
                           [--black_list BLACK_LIST] [--per_read PER_READ]
                           [--hapratio HAPRATIO]
                           [--min_base_quality MIN_BASE_QUALITY]
                           [--mapping_quality MAPPING_QUALITY]
                           [--min_snv MIN_SNV]
                           [--min_read_number MIN_READ_NUMBER]
                           [--min_cg MIN_CG]
                           [--meth_difference METH_DIFFERENCE]
                           [--cpg_difference CPG_DIFFERENCE]
                           [--methyl_coverage METHYL_COVERAGE]
                           [--threads THREADS] [--chunk_size CHUNK_SIZE]
                           [--include_supplementary]

Phasing reads and Methylation

required arguments:
  --bam BAM, -b BAM     The path to the cordinate sorted bam file.
  --output OUTPUT, -o OUTPUT
                        The path to directory and prefix to save files. e.g
                        path/to/directory/prefix
  --vcf VCF, -v VCF     The path to the vcf file.
  --strand_vcf PHASED_VCF, -sv STRAND_VCF
                        The path to the chromosome-scale strand-seq phased vcf file.
 
required arguments if PofO needs to be determined.:
  --known_dmr KNOWN_DMR, -kd KNOWN_DMR
                        The path to the input file for known imprinted
                        DMRs.File must have the following information the
                        following column order: chromosome start end
                        MethylatedAlleleOrigin where origine is the methylated
                        allele origine which must be either maternal or
                        paternal. By default, we use version 1 list in repo's
                        patmat directory.
  --methylcallfile METHYLCALLFILE, -mc METHYLCALLFILE
                        If you want to phase methyl call file (methycall
                        output format) to also calculate methylation frequency
                        for each haplotype give the path to the bgziped
                        methylation call file from methyl_call_processor
                        Module.

Optional arguments.:
  -h, --help            show this help message and exit
  --whatshap_vcf WHATSHAP_VCF, -wv WHATSHAP_VCF
                        Path to the WhatsHap phased vcf file that is produced
                        from phasing nanopore reads using WhatsHap. This can
                        be useful when the chromosome-scale phased variants
                        are very sparce. File must be sorted and indexed using
                        tabix.
  --whatshap_block WHATSHAP_BLOCK, -wb WHATSHAP_BLOCK
                        Path to the WhatsHap block file file. This file can
                        becreated using whatshap stats command. File must
                        beconverted to a bed format with chromosome start end
                        in the first three columns. If no block file is given
                        then the assumption is that the last part after : sign
                        in the 10th column is the phase set (PS) name and
                        blocks will be calculated internaly.
  --black_list BLACK_LIST, -bl BLACK_LIST
                        List of regions to ignore ther strand-seq
                        phasedcstatus three first columns must be chromosome
                        start end. If black list is given the vcf file must be
                        indexed using tabix.
  --per_read PER_READ, -pr PER_READ
                        If it is your second try and you have per read info
                        file from the first try there is no need to give vcf
                        file, instead give the path to the per read info file.
                        This will be significantly faster.
  --hapratio HAPRATIO, -hr HAPRATIO
                        0-1 . Minimmum ratio of variants a read must have from
                        a haplotype to assign it to that haplotype. Default is
                        0.75.
  --min_base_quality MIN_BASE_QUALITY, -mbq MIN_BASE_QUALITY
                        Only include bases with phred score higher or equal
                        than this option. Default is >=7.
  --mapping_quality MAPPING_QUALITY, -mq MAPPING_QUALITY
                        An integer value to specify thereshold for filtering
                        reads based om mapping quality. Default is >=20
  --min_snv MIN_SNV, -ms MIN_SNV
                        minimum number of phased SNVs must a read have to be
                        phased. Default= 1
  --min_read_number MIN_READ_NUMBER, -mr MIN_READ_NUMBER
                        minimum number of reads to support a variant to assign
                        to each haplotype. Default= 2
  --min_cg MIN_CG, -mcg MIN_CG
                        Minimmum number of CpGs an iDMR must have to consider
                        it for PofO assignment. Default is 11.
  --meth_difference METH_DIFFERENCE, -md METH_DIFFERENCE
                        Methylation difference cutoff for HP1-HP2 or HP2-HP1
                        CpG methylation. Default is 0.35.
  --cpg_difference CPG_DIFFERENCE, -cd CPG_DIFFERENCE
                        Cut off for the fraction of CpGs between haplotypes
                        must be differentially methylated at an iDMR to
                        consider it for PofO assignment. Default is 0.1.
  --methyl_coverage METHYL_COVERAGE, -mcov METHYL_COVERAGE
                        Minimmum Coverage at each CpG site when calculating
                        methylation frequency. Default is 1.
  --threads THREADS, -t THREADS
                        Number of parallel processes. Default is 4.
  --chunk_size CHUNK_SIZE, -cs CHUNK_SIZE
                        Chunk per process. Default is 100
  --include_supplementary, -is
                        Also include supplementary reads.
```
