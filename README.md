# PatMat:  
![](docs/FlowChart.png)  

This workflow enables simultanous chromosome-scale haplotyping and parent-of-origin detection in a single sample without any parental data using a combination of nanopore 
sequencing and Strand-seq.
We will use nanopore-detected variants and their long-range phasing from Strand-seq to detect chromosome-scale haplotypes. We then use DNA methylation information at known 
imprinted regions to detect parent-of-origin.  

**Citation:** [Parent-of-origin detection and chromosome-scale haplotyping using long-read DNA methylation sequencing and 
Strand-seq](https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00191-4)  
  
To run this workflow you will need the following third-party tools for processing nanopore data:  
**Guppy**: For basecalling nanopore reads.  
**[Minimap2](https://github.com/lh3/minimap2)**: To align nanopore reads to the reference genome.  
**[Nanopolish](https://github.com/jts/nanopolish)**: To call DNA methylation from nanopore data.  
**[NanoMethPhase](https://github.com/vahidAK/NanoMethPhase)**: To process methylation call results from nanopolish.  
**[Clair3](https://github.com/HKU-BAL/Clair3)**: To call variants from aligned nanopore reads.  

Additional software is required to use Strand-seq, as described [below](https://github.com/vahidAK/PatMat/blob/main/README.md#2--Strand\-seq-Data-Analysis).

The workflow was developed using the above tools. However, you may use alternatives to each tool.    

Finally you need to use our tool in this repository "**PatMat.py**" to detect chromosome-scale parent-of-origin resolved haplotypes. Try to use the latest released version. 
You can download the latest release and unzip or untar the files and use the PatMat.py script in the patmat folder. Alternatively, you can clone the GitHub repository and 
use the PatMat.py in the patmat folder. Before using PatMat.py you need to satisfy (install) the following dependencies:  
[bgzip](http://www.htslib.org/doc/bgzip.html)  
[tabix](http://www.htslib.org/doc/tabix.html)  
python>=3.7.4 and its following dependencies:  
pytabix>=0.1  
pysam>=0.16.0  
tqdm>=4.54.1   

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
Here we use [minimap2](https://github.com/lh3/minimap2) to align nanopore reads to reference genome. You may use other tools such as 
[Winnowmap](https://github.com/marbl/Winnowmap).
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

### 1-3 Indexing and Methylation Calling from nanopore data  
Here we use [nanopolish](https://github.com/jts/nanopolish) for methylation calling however you may use [f5c](https://github.com/hasindu2008/f5c), 
[megalodon](https://github.com/nanoporetech/megalodon) or [DeepSignal](https://github.com/bioinfomaticsCSU/deepsignal). 

#### 1-3-1 indexing fastq file using fast5 files:

NOTE: Fastqs must be merged to a single file

```
nanopolish index -d /path/to/Fast5_files reads.fastq
```
You can also specify sequencing summary file to accelerate indexing.  
Nanopolish index proccess can be time consuming. [f5c](https://github.com/hasindu2008/f5c) which is an optimised and GPU accelerated version of nanopolish can be used for 
indexing fastq using fast5 files on multiple threads. The output of f5c index (for fast5) is [equivalent](https://hasindu2008.github.io/f5c/docs/commands) to that from 
nanopolish index.  

```
f5c index -t <# of threads> --iop <# of I/O processes to read fast5 files> -d /path/to/Fast5_files reads.fastq
```

#### 1-3-2 Methylation calling for CpG from each read:

```
nanopolish call-methylation \
  -t <number_of_threads> -q cpg \
  -r /path/to/Nanopore_reads.fastq \
  -b /path/to/Nanopore_aligned_reads.bam \
  -g /path/to/reference.fa > /path/to/MethylationCall.tsv
```
f5c (versions >=v0.7) can be also used for methylation calling. f5c versions >=v0.7 outputs similar columns as later nanopolish versions (as follows), therefore it is 
compatible with NanoMethPhase.  

```
chromosome	strand	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_motifs	sequence
```
  
#### 1-3-3 Pre-processing methylation call file
We then need to pre-process methylation call file from nanopolish using [NanoMethPhase](https://github.com/vahidAK/NanoMethPhase) methyl_call_processor module.
```
nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > NanoMethPhase_MethylationCall.bed.gz && tabix -p bed 
NanoMethPhase_MethylationCall.bed.gz
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
### 2-1 Quick overview
Typically, 30-100 good quality Strand-seq libraries with at least 20 million unique reads in total should be used for phasing. These must be aligned to the GRCh38 reference 
genome and poor-quality libraries must be identified and removed using [ASHLEYS QC](https://github.com/friendsofstrandseq/ashleys-qc) (Gros et al. 2021).

The inversion-aware Strand-seq phasing routine described here has the following dependencies. Some of these can be installed with the [conda environment 
file](https://github.com/vahidAK/PatMat/tree/main/Strand-seq/env.yml).

* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [R](https://www.r-project.org/) (v4.3.0 or higher)
* The R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
* The R package [BiocManager](https://cran.r-project.org/web/packages/BiocManager/index.html)
* The R package [argparse](https://cran.r-project.org/web/packages/argparse/index.html)
* The R package [BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html)
* The R package [InvertypeR](https://github.com/vincent-hanlon/InvertypeR)

With this software installed and with the BAM files and a nanopore-derived VCF file of non-phased SNVs, something like the following can be run to perform Strand-seq 
phasing (with a Linux OS): 

```
./Strand-seq/strandseq_phase.R \
    -p TRUE \
    -i /path/to/strandseq/bams/ \
    -o ./phased \
    -t 12 \
    -n HG005 \
    /path/to/VCF/of/snvs.vcf
```

(Note that as of June 2023, the StrandPhaseR dependency issues several warning messages ("closing unused connection") when this is run, but that seems to be a bug in the 
dependency rather than an issue with strandseq_phase.R)

It may be necessary to change the permissions first:

```
chmod 770 ./Strand-seq/strandseq_phase.R
```

This method calls inversions using [InvertypeR](https://github.com/vincent-hanlon/InvertypeR) (most of the runtime), which help refine phasing, and then it phases the SNVs 
using the standard Strand-seq R packages [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) and 
[StrandPhaseR](https://github.com/daewoooo/StrandPhaseR). The result is a phased VCF file of SNVs ("samplename.phased.inv_aware.vcf"), which can be used with PatMat.py as 
described below.

Here is the full list of options for `strandseq_phase.R`:

```
usage: ./strandseq_phase.R [-h] [-p TRUE or FALSE] [-i /path/to/input/]
                           [-o /path/to/output/] [-t integer] [-n string]
                           [--inversion_list /path/to/BED]
                           [--hard_mask /path/to/BED]
                           [--soft_mask /path/to/BED]
                           [--prior comma-separated numbers]
                           [--chromosomes comma-separated strings]
                           /path/to/snps.vcf

Performs inversion-aware Strand-seq phasing of a VCF file of SNVs. Requires
bcftools (samtools.github.io/bcftools/bcftools.html), R>=4.3.0, and the R
packages devtools (CRAN), BiocManager (CRAN), InvertypeR (GitHub: vincent-
hanlon/InvertypeR), argparse (CRAN), and BSgenome.Hsapiens.UCSC.hg38
(Bioconductor).

positional arguments:
  /path/to/snps.vcf     Absolute path to a VCF file of SNVs to phase.

options:
  -h, --help            show this help message and exit
  -p TRUE or FALSE, --paired TRUE or FALSE
                        Are the Strand-seq reads paired end? Default: TRUE.
  -i /path/to/input/, --input_folder /path/to/input/
                        Absolute path to the directory containing good-quality
                        Strand-seq libraries for your sample. Default: '.'.
  -o /path/to/output/, --output_folder /path/to/output/
                        Absolute path to a directory where output files should
                        be written. Default: '.'.
  -t integer, --threads integer
                        The number of parallel threads to run on. Default: 4.
  -n string, --name string
                        The name of the sample. This will appear in the VCF
                        file. E.g., 'HG005'. Default: 'unknown'.
  --inversion_list /path/to/BED
                        Absolute path to a BED file containing genomic
                        intervals that might be inversions. This is typically
                        a list from the literature so the file
                        hanlon_2021_BMCgenomics_augmented.bed on the PatMat
                        GitHub (originally from vincent-hanlon/InvertypeR) is
                        a good start. Default: the file suggested above, if it
                        is in the same directory as strandseq_phase.R.
  --hard_mask /path/to/BED
                        Absolute path to a BED file containing regions with
                        unreliable Strand-seq data. The file
                        hard_mask.GRCh38.humans.bed on the PatMat GitHub is a
                        good start. Default: the file suggested above, if it
                        is in the same directory as strandseq_phase.R.
  --soft_mask /path/to/BED
                        Absolute path to a BED file containing regions, like
                        very large inversions, that occasionally interfere
                        with composite file creation. Rarely really necessary
                        (see InvertypeR documentation) Default: the file
                        suggested above, if it is in the same directory as
                        strandseq_phase.R, which contains the three largest
                        autosomal inversions according to Porubsky et al.
                        2022, except for a very rare one on chr2.
  --prior comma-separated numbers
                        A comma-separated list of prior weights for inversion
                        genotypes, without spaces. Only needs to be altered if
                        --inversion_list is not the default. E.g.,
                        0.96,0.02,0.02. See InvertypeR for more details.
  --chromosomes comma-separated strings
                        A comma-separated list of chromosome names, without
                        spaces. Only SNVs on these chromosomes will be phased.
                        E.g., chr1,chr2,chr3. Default: the 22 autosomes.
```
### 2-2 Installations and setup
Note that for the [parent-of-origin phasing paper](https://doi.org/10.1016/j.xgen.2022.100233) that first presented this method, we used the dependencies and code from PatMat v1.1.1 (this repo). However, the most recent version is preferred. 

First, R (v4.3.0 or higher) and bcftools can be installed separately, or using miniconda3 and the "env.yml" file:

```
conda env create --file ./Strand-seq/env.yml -n strandseq
conda activate strandseq
```
The R packages devtools, BiocManager, InvertypeR, argparse, and BSgenome.Hsapiens.UCSC.hg38 must be installed using `install.packages()` from base 
R, `BiocManager::install()`, or `devtools::install_github()`. 

### 2-2-1 Library QC
Separately, the Strand-seq library QC tool [ASHLEYS QC](https://github.com/friendsofstrandseq/ashleys-qc) should be used to select only good-quality libraries for analysis 
(installation instructions in the GitHub link). Typically, after activating the correct conda environment (`conda activate ashleys`), move to the directory containing 
aligned and indexed Strand-seq BAM files and run something like the following:
```
ashleys.py -j 12 features -f ./ -w 5000000 2000000 1000000 800000 600000 400000 200000 -o ./features.tsv
ashleys.py predict -p ./features.tsv -o ./quality.txt -m scripts/tools/svc_default.pkl
```
Then examine quality.txt and either (i) use libraries with a score >0.5 or (ii) use libraries with a score >0.7 and have a domain expert manually inspect libraries with a 
score 0.5-0.7 to see whether they should be included in the analysis. For instructions on how to align FASTQ files, mark duplicates etc., and generate BAM files, see 
alignment.sh in Methods part 1 of this [book chapter](https://dx.doi.org/10.14288/1.0406302), also found [here](https://github.com/vincent-hanlon/MiMB-StrandPhaseR) 
(otherwise, just use your standard sequence alignment for short read data, keeping the single-cell libraries separate).

### 2-3-1 Composite files for inversion calling

In practice, performing Strand-seq phasing now just requires running `./Strand-seq/strandseq_phase.R` as described above. The rest of Section 2 of this user guide just 
describes what the phasing method is actually doing.

Strand-seq libraries are typically low-coverage, and this makes it hard to discover and genotype small inversions. To address this, we use the R package InvertypeR to 
combine data from many libraries into two composite files: one built from regions of libraries where all reads mapped with the same orientation (Watson-Watson or 
Crick-Crick regions), and one where reads mapped with both orientations (Watson-Crick). The latter file entails a phasing step to distinguish cases where (for an autosome) 
homolog 1 gave the forward reads and homolog 2 gave the reverse reads, rather than homolog 1 giving reverse reads and homolog 2 giving forward reads. Apart from identifying 
and phasing such regions, this process is effectively a problem of merging BAM files (loaded into R) and reorienting reads in some cases. 

The main InvertypeR function, which genotypes inversions, takes as input a list of positions to examine. To obtain a list of putative inversions de novo, we run BreakpointR 
on the composite files three times with different bin sizes and extract the coordinates of short segments of the genome with unexpected read orientations (such as inversions
might give). We combine this with a list of inversions from the literature.

### 2-3-2 Inversion genotyping
InvertypeR counts reads inside the putative inversions by orientation in each composite file. Using a simple Bayesian model of read counts, posterior probabilites for 
inversion genotypes are calculated (effectively a comparison of the actual read count pattern with expected read count patterns for the various genotype and error signals). 
InvertypeR also resizes inversions if the coordinates do not perfectly match the region with reversed read orientations. We take inversions with a posterior probability 
above 95% for either the heterozygous or homozygous genotype. Since InvertypeR is run separately (with different prior probabilities) for the putative inversions obtained 
from the literature or from BreakpointR (above), we then combine them by merging overlapping inversions subject to some constraints. Only inversions larger than 10 kb are 
used to correct phasing.

### 2-3-3 Phasing
All the above inversion calling is used to correct or refine SNV phasing inside inversions (where Strand-seq would otherwise make mistakes). This is important for the iDMRs 
that fall inside inversions, and when variants of interest fall inside inversions. We use StrandPhaseR to phase the nanopore-derived SNVs, and then we correct the phasing 
within inversions using the StrandPhaseR tool `correctInvertedRegionPhasing()`. For this process, first we identify Watson-Crick regions (aka WC regions; used for phasing) 
with BreakpointR, and then StrandPhaseR assigns alleles that appear in reads with opposite orientations to different homologs (assuming the reads are in the same WC region 
or chromosome in the same cell). It then combines the phase information from many cells to produce a consensus phased VCF. The inversion correction step then effectively 
switches the haplotypes of alleles inside homozygous inversions and re-phases alleles inside heterozygous inversions. A more complete step-by-step guide to using 
StrandPhaseR (excluding the inversion correction) can be found [here](https://dx.doi.org/10.14288/1.0406302).

## 3- Parent-of-origin detection
Finally, parent-of-origin chromosome-scale haplotypes can be built using PatMat.py:  
```
PatMat.py -v /path/to/Passed_Clair3_Variants.vcf \
 -sv /path/to/StrandSeq_phased_variants.vcf \
 -mc /path/to/NanoMethPhase_MethylationCall.bed.gz \
 -b /path/to/Nanopore_aligned_reads.bam \
 -o <Output pass and prefix> \
 -t <# of threads>
```

Here is the full list of options:  
```
python patmat/PatMat.py -h

usage: PatMat.py [-h] --bam BAM --output OUTPUT --vcf VCF --strand_vcf
                 STRAND_VCF --methylcallfile METHYLCALLFILE
                 [--known_dmr KNOWN_DMR] [--whatshap_vcf WHATSHAP_VCF]
                 [--whatshap_block WHATSHAP_BLOCK] [--black_list BLACK_LIST]
                 [--per_read PER_READ] [--hapratio HAPRATIO]
                 [--min_base_quality MIN_BASE_QUALITY]
                 [--mapping_quality MAPPING_QUALITY]
                 [--min_variant MIN_VARIANT]
                 [--min_read_number MIN_READ_NUMBER] [--min_cg MIN_CG]
                 [--meth_difference METH_DIFFERENCE]
                 [--cpg_difference CPG_DIFFERENCE]
                 [--methyl_coverage METHYL_COVERAGE] [--threads THREADS]
                 [--chunk_size CHUNK_SIZE] [--include_supplementary]
                 [--include_indels] [--version]

Phasing reads and Methylation using strand-seq and nanopore to determine PofO
of each homologous chromosome in a single sample.

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  --bam BAM, -b BAM     The path to the cordinate sorted bam file.
  --output OUTPUT, -o OUTPUT
                        The path to directory and prefix to save files. e.g
                        path/to/directory/prefix
  --vcf VCF, -v VCF     The path to the vcf file.
  --strand_vcf STRAND_VCF, -sv STRAND_VCF
                        The path to the chromosome-scale phased vcf file.This
                        is the input vcf file that has been phased using
                        strand-seq data.
  --methylcallfile METHYLCALLFILE, -mc METHYLCALLFILE
                        The path to the bgziped and indexed methylation call
                        file processed using NanoMethPhase
                        methyl_call_processor Module.

Optional arguments.:
  --known_dmr KNOWN_DMR, -kd KNOWN_DMR
                        The path to the input file for known imprinted DMRs.
                        File must have the following information in the
                        following column order: chromosome start end
                        MethylatedAlleleOrigin where the methylated allele
                        origin must be either maternal or paternal (First row
                        must be header). By default, we use iDMR list in
                        repo's patmat directory.
  --whatshap_vcf WHATSHAP_VCF, -wv WHATSHAP_VCF
                        Path to the WhatsHap phased vcf file that is produced
                        from phasing input vcf file using nanopore reads via
                        WhatsHap. This can be useful when the chromosome-scale
                        phased variants are very sparce. File must be sorted
                        and indexed using tabix.
  --whatshap_block WHATSHAP_BLOCK, -wb WHATSHAP_BLOCK
                        Path to the WhatsHap block file. This file can be
                        created using whatshap stats command. File must be
                        converted to a bed format with chromosome start end in
                        the first three columns (First row must be header). If
                        no block file is given then the assumption is that the
                        last part after : sign in the 10th column is the phase
                        set (PS) name and blocks will be calculated
                        internally.
  --black_list BLACK_LIST, -bl BLACK_LIST
                        List of regions to ignore phased varinats at them.
                        Three first columns must be chromosome start end. If
                        black list is given the vcf file must be indexed using
                        tabix.
  --per_read PER_READ, -pr PER_READ
                        If it is your second try and you have per read info
                        file give the path to the per read info file. This
                        will be significantly faster. This is also useful when
                        you want to try different thresholds for options (Note
                        that if you also provided WhatsHap phased vcf in your
                        first try, then you cannot use per-read to try
                        different --min_variant or --hapratio because these
                        options will be also used to correct WhatsHap phased-
                        block switches using strand-seq phased variants),
                        different dmr list, black list, include/exclude
                        indels, and include/exclude supp reads.
  --hapratio HAPRATIO, -hr HAPRATIO
                        0-1 . Minimum ratio of variants a read must have from
                        a haplotype to assign it to that haplotype. Default is
                        0.75. Note that if you also provide WhatsHap phased
                        vcf file this option will be also used to correct
                        phased-block switches using Strand-seq phased
                        variants. In this case, it is minimum ratio of phased
                        variants at a block that supports the dicision based
                        on strand-seq phased varinats.
  --min_base_quality MIN_BASE_QUALITY, -mbq MIN_BASE_QUALITY
                        Only include bases with phred score higher or equal
                        than this option. Default is >=7.
  --mapping_quality MAPPING_QUALITY, -mq MAPPING_QUALITY
                        An integer value to specify thereshold for filtering
                        reads based on mapping quality. Default is >=20
  --min_variant MIN_VARIANT, -mv MIN_VARIANT
                        Minimum number of phased variants must a read have to
                        be phased. Default= 1. Note that if you also provide
                        WhatsHap phased vcf file this option will be also used
                        to correct phased-block switches using Strand-seq
                        phased variants. In this case, it is the minimum
                        number of phased variants at a block that need to
                        support the dicision based on strand-seq phased
                        varinats.
  --min_read_number MIN_READ_NUMBER, -mr MIN_READ_NUMBER
                        Minimum number of reads to support a variant to assign
                        to each haplotype. Default= 2
  --min_cg MIN_CG, -mcg MIN_CG
                        Minimum number of CpGs an iDMR must have to consider
                        it for PofO assignment. Default is 12.
  --meth_difference METH_DIFFERENCE, -md METH_DIFFERENCE
                        0-1. Minimum methylation difference cutoff for HP1-HP2
                        or HP2-HP1 CpG methylation. Default is 0.35.
  --cpg_difference CPG_DIFFERENCE, -cd CPG_DIFFERENCE
                        0-1. Minimum cut off for the fraction of CpGs between
                        haplotypes must be differentially methylated at an
                        iDMR to consider it for PofO assignment. Default is
                        0.1.
  --methyl_coverage METHYL_COVERAGE, -mcov METHYL_COVERAGE
                        Minimum Coverage at each CpG site when calculating
                        methylation frequency. Default is 1.
  --threads THREADS, -t THREADS
                        Number of parallel processes. Default is 4.
  --chunk_size CHUNK_SIZE, -cs CHUNK_SIZE
                        Chunk per process. Default is 100
  --include_supplementary, -is
                        Also include supplementary reads (Not recommended).
  --include_indels, -ind
                        Also include indels for read phasing to haplotypes.
  --version             show program's version number and exit
```
### 3-1- Outputs
PatMat will generate multiple outputs.
#### 3-1-1 NonPofO_HP1-HP2 (Non parent-of-origin) results 
These are a vcf and a tsv file. These files represent the results for phasing reads and re-phasing het variants (haplotype 1 or HP1 and haplotype 2 or HP2) before assigning 
parent-of-origin. In the vcf file, for phased 0/1 (0|1 or 1|0) variants the last column includes HP1|HP2 (Ref is HP1 and alt is HP2), or HP2|HP1 (Ref is HP2 and alt is HP1) 
and for the phased 1/2 variants (1|2) the last column includes Ref_HP1|HP2 (the part before comma on the 5th column is HP1 and the part after comma is HP2) or Ref_HP2|HP1 
(the part before comma on the 5th column is HP2 and the part after comma is HP1).  
Note: During re-phasing input vcf variants, if your input vcf is a phased vcf file and some of the varinats could not be re-phased, the phase sign "|" will be just replaced 
by "/" sign (e.g. 1|0 will be 1/0).  
#### 3-1-2 PofO_Assignment results 
These are a vcf and a tsv file. These files represent the results after assigning the parent-of-origin to HP1 and HP2 reads and variants. In the vcf file, for phased 0/1 
(0|1 or 1|0) variants the last column includes Mat|Pat (Ref is maternal and alt is paternal), or Pat|Mat (Ref is paternal and alt is maternal) and for the phased 1/2 
variants (1|2) the last column includes Ref_Mat|Pat (the part before comma on the 5th column is maternal and the part after comma is paternal) or Ref_Pat|Mat (the part 
before comma on the 5th column is paternal and the part after comma is maternal).  
Note: During PofO assignment to the re-phased variants the phase sign "|" will be just replaced by "/" sign (e.g. 1|0 will be 1/0) if PofO could not be inferred.  
#### 3-1-3 CpG-Methylation-Status-at-DMRs 
This file represents status of CpGs and their methylation at each DMR on each haplotype, including the number or common CpGs between haplotypes and their methylation 
frequencies on each haplotypes, how many of them showed given methylation difference on each haplotype, and contribution or detection value of the DMR for each haplotype.  
#### 3-1-4 HP1_HP2_PerReadInfo 
This file includes per-read information including coordinates and strand of the reads on reference, read IDs, read flag and if the read is supplemenary or not, read mapping 
quality, and finally the positions, phred score base quality and base(s) from the read at the input phased het variants from strand-seq (or strand-seq plus WhatsHap, if 
given) for HP1 and HP2 and unphased variants in the input vcf file (1|2 varinats are considered as unphased). Base quality for indels represent the base quality of the 
first base. Positions in the per-read file are zero-based.  

**Note:** If you wish to try different criteria, the per-read file produced by PatMat >=v1.2.0 allows you to try different thresholds for options (**Note** that if you also 
provided WhatsHap phased vcf in your first try, then you **cannot** use per-read to try different --min_variant or --hapratio because these options will be also used to 
correct WhatsHap phased-block switches using strand-seq phased variants.), different dmr list, black list, include/exclude indels, and include/exclude supp reads much 
faster. These are also true for previous versions and their per-read **except** that per-read from previous versions **cannot** be used for different black list or 
include/exclude supp reads. NanoMethPhase phase module also produces a per-read file, however, per-read file from PatMat is NOT equivalent to the per-read file from 
NanoMethPhase. 
