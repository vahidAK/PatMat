# PatMat:  
![](docs/FlowChart.png)  

This workflow enables simultaneous chromosome-scale haplotyping and parent-of-origin detection in a single sample without any parental data using a combination of nanopore 
sequencing and Strand-seq.
We will use nanopore-detected variants and their long-range phasing from Strand-seq to detect chromosome-scale haplotypes. We then use DNA methylation information at known 
imprinted regions to detect parent-of-origin.  

**Citation:** [Parent-of-origin detection and chromosome-scale haplotyping using long-read DNA methylation sequencing and 
Strand-seq](https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00191-4)  

Table of Contents
=================
* [Installation](https://github.com/vahidAK/PatMat/blob/main/README.md#installation)
* [Full Tutorial](https://github.com/vahidAK/PatMat/blob/main/README.md#full-tutorial)
  * [Nanopore Data Analysis](https://github.com/vahidAK/PatMat/blob/main/README.md#1--nanopore-data-analysis)
    * [Basecalling, mapping, and methylation calling from nanopore data using Guppy](https://github.com/vahidAK/PatMat/blob/main/README.md#1\-1-basecalling-mapping-and-methylation-calling-from-nanopore-data-using-guppy)
    * [Variant Calling from nanopore data using clair3](https://github.com/vahidAK/PatMat/blob/main/README.md#1\-2-variant-calling-from-nanopore-data-using-clair3)
  * [Strand-seq Data Analysis](https://github.com/vahidAK/PatMat/blob/main/README.md#2--strand\-seq-data-analysis)
    * [Library QC](https://github.com/vahidAK/PatMat/blob/main/README.md#2\-1-library-qc)
    * [Phasing](https://github.com/vahidAK/PatMat/blob/main/README.md#2\-2-phasing)
    * [What is actually going on here?](https://github.com/vahidAK/PatMat/blob/main/README.md#2\-3-what-is-actually-going-on-here)
      * [Composite files](https://github.com/vahidAK/PatMat/blob/main/README.md#2\-3-1-composite-files)
      * [Inversion genotyping](https://github.com/vahidAK/PatMat/blob/main/README.md#2\-3-2-inversion-genotyping)
      * [Phasing](https://github.com/vahidAK/PatMat/blob/main/README.md#2\-3-3-phasing)
  * [Parent-of-origin detection](https://github.com/vahidAK/PatMat/blob/main/README.md#3--parent\-of\-origin-detection)
    * [Outputs](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-outputs)
      * [NonPofO_HP1-HP2 results](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-1-nonpofo_hp1\-hp2-results)
      * [PofO_Assignment results](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-2-pofo_assignment-results)
      * [CpG-Methylation-Status-at-DMRs](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-3-cpg\-methylation\-status\-at\-dmrs)
      * [DMLtest.tsv.gz and callDML.tsv.gz](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-4-dmltesttsvgz-and-calldmltsvgz)
      * [PofO_Scores.tsv](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-5-pofo_scorestsv)
      * [HP1_HP2_PerReadInfo](https://github.com/vahidAK/PatMat/blob/main/README.md#3-1-6-hp1_hp2_perreadinfo)
* [More info about other methylation callers](https://github.com/vahidAK/PatMat/blob/main/README.md#more-info-about-other-methylation-callers)
        
  
# Installation
The workflow is basically two parts, the nanopore analysis part and the Strand-seq analysis part. All the tools needed to run the workflow can be installed by downloading/cloning this repository as explained below. However, there are a few dependencies that you need to have/install separately first including conda, guppy, clair3, and ashleys (See the notes below).   
**Note 1**: You first need to have conda/miniconda installed. If you do not have conda/miniconda follow the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install it and then run the following commands. Alternatively, [mamba](https://mamba.readthedocs.io/en/latest/installation.html) is a faster alternative to conda.   
**Note 2**: For the nanopore analysis part you need guppy basecaller that does basecalling, mapping, and methylation calling. This tool is only available through the [Oxford Nanopore Technologies community website](https://nanoporetech.com/community). You need to create an account there and get this tool yourself.  
**Note 3**: For variant calling from nanopore data you need [clair3](https://github.com/HKU-BAL/Clair3). You need to install clair3 in its own dedicated conda environment by running ```conda create -n clair3 -c bioconda clair3 python=3.9.0 -y```. For more instructions and also to download the appropriate model for variant calling visit [clair3 GitHub](https://github.com/HKU-BAL/Clair3).  
**Note 4**: For the quality control of Strand-seq libraries you need [ASHLEYS QC](https://github.com/friendsofstrandseq/ashleys-qc). It needs to be installed in its own dedicated conda environment. Follow their instructions on [ASHLEYS QC GitHub](https://github.com/friendsofstrandseq/ashleys-qc) to install this tool.  

To use the tools in this repository and run the workflow you can download the latest release or clone the repository and install the required dependencies in the [env.yml](https://github.com/vahidAK/PatMat/blob/main/env.yml) as follow:  

**NOTE:** You must have conda and mamba installed before installation.

To clone and install:
```
git clone https://github.com/vahidAK/PatMat.git
cd PatMat
bash -l ./install.sh
# Now we need to open R and install InvertypeR:
conda activate patmat
R
devtools::install_github("vincent-hanlon/InvertypeR")

``` 
OR to download the latest release and install:
```
VERSION=1.3.0
wget https://github.com/vahidAK/PatMat/archive/refs/tags/v"$VERSION".tar.gz && tar -xzf v"$VERSION".tar.gz
cd PatMat-"$VERSION"/
bash -l ./install.sh
# Now we need to open R and install InvertypeR:
conda activate patmat
R
devtools::install_github("vincent-hanlon/InvertypeR")

```
The above commands will clone/download the repository and install all the dependencies in the patmat environment. You need to first activate the environment to be able to run the tools  ```conda activate patmat```.

The `pak` installations require that you have valid GitHub credentials on your system or no GitHub credentials.

# Full Tutorial  
Note that for the [parent-of-origin phasing paper](https://doi.org/10.1016/j.xgen.2022.100233) that first presented this method, we used the dependencies and code from PatMat v1.1.1 (this repo). However, the most recent version is preferred.  
Currently, the workflow is available for the human reference genome GRCh38/hg38 because iDMR coordinates are based on hg38 and Strand-seq analysis is also configured based on hg38.   

## 1- Nanopore Data Analysis
### 1-1 Basecalling, mapping, and methylation calling from nanopore data using Guppy

```
guppy_basecaller  \
    --input_path <Path to fast5 directory> \
    --save_path <path to output directory> \
    --config <appropriate configuration file (e.g. dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_sup_prom.cfg)> \
    --device <GPU devices to be used (e,g. cuda:0 cuda:1)> \
    --compress_fastq \
    --recursive \
    --bam_out \
    --index \
    --align_ref ${reference_fasta} \
    --disable_pings \
    --data_path <Path to use for loading any data files the application requires (e.g., config file). This is needed when your given config file is not among the default configs shipped by the tool> \ 
    --trim_strategy dna
```
After basecalling if you have multiple bam files you need to merge them all into a single bam file. Bam file must be reference coordinate sorted and indexed.  

### 1-2 Variant calling from nanopore data using clair3
Activate the conda environment:
```
conda activate clair3
```
Now you can call variants using clair3:
```
run_clair3.sh --bam_fn=/path/to/Nanopore_aligned_reads.bam \
  --ref_fn=/path/to/reference.fa \
  --output=/path/to/output/directory \
  --threads=<# of threads> --platform=ont \
  --model_path=/path/to/model/ont_guppy5_r941_sup_g5014
```
After variant calling the results will be in the merge_output.vcf.gz file in the output directory.  

## 2- Strand-seq Data Analysis

Typically, 30-100 good quality Strand-seq libraries with at least 20 million unique reads in total should be used for phasing. These must be aligned to the GRCh38 reference
genome and poor-quality libraries must be identified and removed using [ASHLEYS QC](https://github.com/friendsofstrandseq/ashleys-qc) (Gros et al. 2021). Then 
`strandseq_phase.R` can be run from the command line to call inversions and perform inversion-aware phasing.

### 2-1 Library QC
The Strand-seq library QC tool [ASHLEYS QC](https://github.com/friendsofstrandseq/ashleys-qc) should be used to select only good-quality libraries for analysis
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

### 2-2 Phasing

Activate the conda environment:
```
conda activate patmat
```
Then, with the Strand-seq BAM files and a nanopore-derived VCF file of non-phased SNVs, something like the following can be run to perform Strand-seq 
phasing (from the Linux command line):  

```
./Strand-seq/strandseq_phase.R \
    -p TRUE \
    -i /path/to/strandseq/bams/ \
    -o ./phased \
    -t 12 \
    -n samplename \
    /path/to/VCF/of/snvs.vcf
```
Note that as of June 2023, the StrandPhaseR dependency sometimes issues several warning messages ("closing unused connection") when this is run, but that seems to be a bug 
in the dependency rather than an issue with `strandseq_phase.R`. 

The result is a phased VCF file of SNVs ("samplename.phased.inv_aware.vcf"), which can be used with `patmat.py` as 
described below.

Here is the full list of options for `strandseq_phase.R`:

```
usage: ./strandseq_phase.R [-h] [-p TRUE or FALSE] [-i /path/to/BAMs/]
                           [-o /path/to/output/] [-t integer] [-n string]
                           [--inversion_list /path/to/BED]
                           [--hard_mask /path/to/BED]
                           [--soft_mask /path/to/BED]
                           [--prior comma-separated numbers]
                           [--chromosomes comma-separated strings]
                           [--filter_snvs TRUE or FALSE]
                           [--fix_het_invs TRUE or FALSE]
                           /path/to/snvs.vcf

Performs inversion-aware Strand-seq phasing of a VCF file of SNVs. Requires
bcftools (samtools.github.io/bcftools/bcftools.html), R>=4.3.0, and the R
packages InvertypeR (GitHub: vincent-hanlon/InvertypeR), argparse (CRAN), and
BSgenome.Hsapiens.UCSC.hg38 (Bioconductor). These are best installed with the
R package installer, pak (CRAN)!

positional arguments:
  /path/to/snvs.vcf     Absolute path to a VCF file of SNVs to phase.

options:
  -h, --help            show this help message and exit
  -p TRUE or FALSE, --paired TRUE or FALSE
                        Are the Strand-seq reads paired end? Default: TRUE.
  -i /path/to/BAMs/, --input_folder /path/to/BAMs/
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
  --filter_snvs TRUE or FALSE
                        Should SNVs be filtered such that the FILTER column of
                        the input VCF is either 'PASS' or '.', removing
                        potential low quality sites? Default TRUE.
  --fix_het_invs TRUE or FALSE
                        If TRUE, an attempt will be made to correct Strand-seq
                        phasing errors caused by heterozygous inversions. If
                        FALSE, SNVs inside heterozygous inversions will simply
                        not be phased. Switch errors sometimes remain after
                        phase correction at heterozygous inversions, but more
                        SNVs can generally be phased. Default FALSE.
```


### 2-3 What is actually going on here?
In practice, performing Strand-seq phasing now just requires running `./Strand-seq/strandseq_phase.R` as described above. The rest of Section 2 of this user guide just 
describes what the phasing method is actually doing.

#### 2-3-1 Composite files
Strand-seq libraries are typically low-coverage, and this makes it hard to discover and genotype small inversions. To address this, we use the R package InvertypeR to 
combine data from many libraries into two composite files: one built from regions of libraries where all reads mapped with the same orientation (Watson-Watson or 
Crick-Crick regions), and one where reads mapped with both orientations (Watson-Crick). The latter file entails a phasing step to distinguish cases where (for an autosome) 
homolog 1 gave the forward reads and homolog 2 gave the reverse reads, rather than homolog 1 giving reverse reads and homolog 2 giving forward reads. Apart from identifying 
and phasing such regions, this process is effectively a problem of merging BAM files (loaded into R) and reorienting reads in some cases.  

#### 2-3-2 Inversion genotyping
The main InvertypeR function, which genotypes inversions, takes as input a list of positions to examine. To obtain a list of putative inversions de novo, we run BreakpointR 
on the composite files three times with different bin sizes and extract the coordinates of short segments of the genome with unexpected read orientations (such as inversions might give). We combine this with a list of inversions from the literature.

InvertypeR counts reads inside the putative inversions by orientation in each composite file. Using a simple Bayesian model of read counts, posterior probabilites for 
inversion genotypes are calculated (effectively a comparison of the actual read count pattern with expected read count patterns for the various genotype and error signals). 
InvertypeR also resizes inversions if the coordinates do not perfectly match the region with reversed read orientations. We take inversions with a posterior probability 
above 95% for either the heterozygous or homozygous genotype. Since InvertypeR is run separately (with different prior probabilities) for the putative inversions obtained 
from the literature or from BreakpointR (above), we then combine them by merging overlapping inversions subject to some constraints. Only inversions larger than 10 kb are 
used to correct phasing.

#### 2-3-3 Phasing
All the above inversion calling is used to correct or refine SNV phasing inside inversions (where Strand-seq would otherwise make mistakes). This is important for the iDMRs 
that fall inside inversions, and when variants of interest fall inside inversions. We use StrandPhaseR to phase the nanopore-derived SNVs, and then we correct the phasing 
within inversions using the StrandPhaseR tool `correctInvertedRegionPhasing()`. For this process, first we identify Watson-Crick regions (aka WC regions; used for phasing) 
with BreakpointR, and then StrandPhaseR assigns alleles that appear in reads with opposite orientations to different homologs (assuming the reads are in the same WC region 
or chromosome in the same cell). It then combines the phase information from many cells to produce a consensus phased VCF. The inversion correction step then effectively 
switches the haplotypes of alleles inside homozygous inversions and either (i) re-phases alleles inside heterozygous inversions or (ii) removes snvs inside heterozygous 
inversions from the output VCF file. A more complete step-by-step guide to using StrandPhaseR (excluding the inversion correction) can be found [here](https://dx.doi.org/10.14288/1.0406302).

## 3- Parent-of-origin detection
Finally, parent-of-origin assigned chromosome-scale haplotypes can be built using patmat.py.  
Activate the conda environment:
```
conda activate patmat
```
Now you can run patmat.py:
```
patmat.py -v /path/to/Passed_Clair3_Variants.vcf \
 -sv /path/to/StrandSeq_phased_variants.vcf \
 -mc /path/to/guppy.bam \
 -b /path/to/guppy.bam \
 -ref /path/to/reference.fa \
 -o <Output pass and prefix> \
 -t <# of threads>
```

Here is the full list of options:  
```
patmat.py -h

usage: patmat.py --bam BAM --output OUTPUT --vcf VCF --strand_vcf STRAND_VCF
                 --methylcallfile METHYLCALLFILE
                 [--tool_and_callthresh TOOL_AND_CALLTHRESH]
                 [--reference REFERENCE] [--known_dmr KNOWN_DMR]
                 [--whatshap_vcf WHATSHAP_VCF]
                 [--whatshap_block WHATSHAP_BLOCK] [--black_list BLACK_LIST]
                 [--hapratio HAPRATIO] [--min_base_quality MIN_BASE_QUALITY]
                 [--mapping_quality MAPPING_QUALITY]
                 [--min_variant MIN_VARIANT]
                 [--min_read_number MIN_READ_NUMBER] [--min_cg MIN_CG]
                 [--cpg_difference CPG_DIFFERENCE] [--include_all_variants]
                 [--include_supplementary] [--include_indels]
                 [--per_read PER_READ] [--processes PROCESSES]
                 [--chunk_size CHUNK_SIZE] [--delta_cutoff DELTA_CUTOFF]
                 [--pvalue PVALUE] [--smoothing_span SMOOTHING_SPAN]
                 [--smoothing_flag SMOOTHING_FLAG] [--equal_disp EQUAL_DISP]
                 [--dss_processes DSS_PROCESSES] [--version] [-h]

Phasing reads and Methylation using strand-seq and nanopore to determine PofO
of each homologous chromosome in a single sample.

Required arguments:
  --bam BAM, -b BAM     The path to the coordinate sorted bam file.
  --output OUTPUT, -o OUTPUT
                        The path to directory and prefix to save files. e.g
                        path/to/directory/prefix
  --vcf VCF, -v VCF     The path to the vcf file.
  --strand_vcf STRAND_VCF, -sv STRAND_VCF
                        The path to the chromosome-scale phased vcf file. This
                        is the input vcf file that has been phased using
                        strand-seq data.
  --methylcallfile METHYLCALLFILE, -mc METHYLCALLFILE
                        The path to the per-read methylation call file or the
                        bam file with methylation tag. If your input bam file
                        includes methylation tags you just need to specify
                        your input bam file again for this option.
  --reference REFERENCE, -ref REFERENCE
                        If you have given a bam file with methylation tag for
                        the --tool_and_callthresh option, then you must also
                        give the path to the reference file. File must be
                        indexed using samtools faidx.

Optional arguments:
  --tool_and_callthresh TOOL_AND_CALLTHRESH, -tc TOOL_AND_CALLTHRESH
                        Software you have used for methylation calling:Call
                        threshold. Supported files include methbam
                        (Methylation bam format produced by guppy basecaller)
                        and per-read CpG methylation calls from nanoplish (or
                        f5c>=v0.7), megalodon, and deepsignal. For example,
                        nanopolish:1.5 is when methylation calling performed
                        by nanopolish and a CpG with llr >= 1.5 will be
                        considered as methylated and llr <= -1.5 as
                        unmethylated, anything in between will be considered
                        as ambiguous call and ignored. For methbam, megalodon,
                        and deepsignl call threshold will be delta probability
                        (0-1). For example threshold 0.4 means any call >=0.7
                        is methylated and <=0.3 is not and between 0.3-0.7
                        will be ignored. Default is methbam:0.4. If methbam is
                        selected you must also provide path to the reference
                        file using --reference option.
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
                        from phasing input vcf file using nanopore reads. This
                        can be useful when the chromosome-scale phased
                        variants are very sparse but be aware that this might
                        also results in more errors. File must be sorted and
                        indexed using tabix.
  --whatshap_block WHATSHAP_BLOCK, -wb WHATSHAP_BLOCK
                        Path to the WhatsHap block file. This file can be
                        created using whatshap stats command. File must be
                        converted to a bed format with chromosome start end in
                        the first three columns (First row must be header). If
                        no block file is given then the assumption is that the
                        last part after ":" sign in the 10th column is the
                        phase set (PS) name and blocks will be calculated
                        internally.
  --black_list BLACK_LIST, -bl BLACK_LIST
                        List of regions to ignore phased variants at them.
                        Three first columns must be chromosome start end. If a
                        black list is given the vcf file must be indexed using
                        tabix.
  --hapratio HAPRATIO, -hr HAPRATIO
                        0-1. Minimum ratio of variants a read must have from a
                        haplotype to assign it to that haplotype. Default is
                        0.75. Note that if you also provide WhatsHap phased
                        vcf file this option will be also used to correct
                        phased-block switches using Strand-seq phased
                        variants. In this case, it is the minimum ratio of
                        phased variants at a block that supports the decision
                        based on strand-seq phased variants.
  --min_base_quality MIN_BASE_QUALITY, -mbq MIN_BASE_QUALITY
                        Only include bases with phred score higher or equal to
                        this option. The default is >=7. if any read/base in
                        alignment file does not have base quality data or
                        cannot be obtained, this option will be ignored for
                        such reads/bases and all the bases will be used.
  --mapping_quality MAPPING_QUALITY, -mq MAPPING_QUALITY
                        An integer value to specify threshold for filtering
                        reads based on mapping quality. Default is >=20
  --min_variant MIN_VARIANT, -mv MIN_VARIANT
                        Minimum number of phased variants must a read have to
                        be phased. Default= 1. Note that if you also provide
                        WhatsHap phased vcf file this option will be also used
                        to correct phased-block switches using Strand-seq
                        phased variants. In this case, it is the minimum
                        number of phased variants at a block that need to
                        support the decision based on strand-seq phased
                        variants.
  --min_read_number MIN_READ_NUMBER, -mr MIN_READ_NUMBER
                        Minimum number of reads to support a variant to assign
                        to each haplotype. Default= 2
  --min_cg MIN_CG, -mcg MIN_CG
                        Minimum number of CpGs an iDMR must have to consider
                        it for PofO assignment. Default is 5.
  --cpg_difference CPG_DIFFERENCE, -cd CPG_DIFFERENCE
                        Minimum cut off for the fraction of CpGs between
                        haplotypes must be differentially methylated at an
                        iDMR to consider it for PofO assignment. Default is
                        0.1.
  --include_all_variants, -iav
                        By default, only variants that have "PASS" or "." in
                        the FILTER column of the input vcf file will be used
                        during phasing and PofO assignment. Select this flag
                        if you want to use all the variants.
  --include_supplementary, -is
                        Also include supplementary reads.
  --include_indels, -ind
                        Also include indels for read phasing to haplotypes.
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
  --processes PROCESSES, -p PROCESSES
                        Number of parallel processes. Default is 4.
  --chunk_size CHUNK_SIZE, -cs CHUNK_SIZE
                        Chunk per process. Default is 100

Optional arguments. The following options are DSS options for differential methylation analysis to find differentially methylated CpGs between haplotypes:
  --delta_cutoff DELTA_CUTOFF, -dc DELTA_CUTOFF
                        0-1. A threshold for defining differentially
                        methylated loci (DML) or CpGs. In DML testing
                        procedure, hypothesis test that the two groups means
                        are equal is conducted at each CpG site. Here if delta
                        is specified, the function will compute the posterior
                        probability that the difference of the means are
                        greater than delta, and then call DML based on that.
                        Default is 0.075.
  --pvalue PVALUE, -pv PVALUE
                        0-1. When delta is not specified, this is the
                        threshold of p-value for defining DML and loci with
                        p-value less than this threshold will be deemed DMLs.
                        When delta is specified, CpG sites with posterior
                        probability greater than 1-pvalue_threshold are deemed
                        DML. Default is 0.001
  --smoothing_span SMOOTHING_SPAN, -sms SMOOTHING_SPAN
                        The size of smoothing window, in basepairs. Default is
                        500.
  --smoothing_flag SMOOTHING_FLAG, -smf SMOOTHING_FLAG
                        TRUE/FALSE. A flag to indicate whether to apply
                        smoothing in estimating mean methylation levels. For
                        more instructions see the DSS R package guide. Default
                        is TRUE.
  --equal_disp EQUAL_DISP, -ed EQUAL_DISP
                        TRUE/FALSE. A flag to indicate whether the dispersion
                        in two groups are deemed equal or not. For more
                        instructions see the DSS R package guide. Default is
                        FALSE Because there is no biological replicate here,
                        you should specify either equal_disp TRUE or
                        smoothing_flag TRUE. Do not specify both as FALSE.
  --dss_processes DSS_PROCESSES, -dp DSS_PROCESSES
                        Number of parallel processes use for DSS differential
                        methylation analysis. If not given, it will be the
                        same as --processes option. Differential methylation
                        analysis usually requires high memory and if there is
                        not enough memory, specify less number processes using
                        this flag for DSS to allocate available memory for
                        less processes.

Help and version options:
  --version             Print program's version and exit
  -h, --help            Print this help and exit.
```
### 3-1 Outputs
PatMat will generate multiple outputs.
#### 3-1-1 NonPofO_HP1-HP2 results 
These are vcf and tsv files. These files represent the results for phasing methylation and reads and re-phasing het variants (haplotype 1 or HP1 and haplotype 2 or HP2) before assigning parent-of-origin. In the vcf file, for phased 0/1 (0|1 or 1|0) variants the last column includes HP1|HP2 (Ref is HP1 and alt is HP2), or HP2|HP1 (Ref is HP2 and alt is HP1) and for the phased 1/2 variants (1|2) the last column includes Ref_HP1|HP2 (the part before comma on the 5th column is HP1 and the part after comma is HP2) or Ref_HP2|HP1 (the part before comma on the 5th column is HP2 and the part after comma is HP1).   
#### 3-1-2 PofO_Assignment results 
These are vcf and tsv files. These files represent the results after assigning the parent-of-origin to HP1 and HP2 methylation data, reads, and variants. In the vcf file, for phased 0/1 (0|1 or 1|0) variants the last column includes Mat|Pat (Ref is maternal and alt is paternal), or Pat|Mat (Ref is paternal and alt is maternal) and for the phased 1/2 variants (1|2) the last column includes Ref_Mat|Pat (the part before comma on the 5th column is maternal and the part after comma is paternal) or Ref_Pat|Mat (the part before comma on the 5th column is paternal and the part after comma is maternal).  
#### 3-1-3 CpG-Methylation-Status-at-DMRs 
This file represents the status of methylation at iDMRs on each haplotype, including the number of CpGs at haplotypes and their methylation frequencies on each haplotype, how many of them showed differential methylation on each haplotype, and if the iDMR included for PofO assignment and score calculation or not.  
If a chromosome could not be assigned a PofO or the PofO score for a chromosome is low you can manually inspect this file to better understand the methylation status at the iDMRs for the chromosome and manually assign PofO to the chromosome or ignore the chromosome.  
#### 3-1-4 DMLtest.tsv.gz and callDML.tsv.gz
These are the files from DSS statistical analysis for the detection of differentially methylated CpGs. DMLtest stores statistical results for all the CpGs and callDML stores differentially methylated CpGs.
#### 3-1-5 PofO_Scores.tsv 
This file includes the PofO assignment score along with some more information for each chromosome. As of version 1.3.0 PofO assignment scoring is changed and for each chromosome it represents (# of all differentially methylated CGs in chromosome supporting PofO assignment) / (# of all differentially methylated CGs in chromosome).  
#### 3-1-6 HP1_HP2_PerReadInfo 
This file includes per-read information including coordinates and strand of the reads on reference, read flag and if the read is supplementary or not, read alignment length and mapping quality, and finally the positions, phred score base quality and base(s) from the read at the input phased het variants from strand-seq (or strand-seq plus WhatsHap, if given) for HP1 and HP2 and unphased variants in the input vcf file (1|2 variants are considered as unphased). Base quality for indels represents the base quality of the first base. Positions in the per-read file are zero-based.  
**Note:** If you wish to try different criteria, the per-read file produced by PatMat >=v1.2.0 allows you to try different thresholds for options (**Note** that if you also provided WhatsHap phased vcf in your first try, then you **cannot** use per-read to try different --min_variant or --hapratio because these options will be also used to correct WhatsHap phased-block switches using strand-seq phased variants.), different dmr list, black list, include/exclude indels, and include/exclude supp reads much 
faster.  
# More info about other methylation callers
By default, we assume that you called methylation using guppy and your guppy bam also has methylation tag. However, patmat.py also supports other methylation callers including Nanopolish (Or f5c>=v0.7, which is a GPU implementation of nanopolish with the same output format), Megalodon, and DeepSignal. Here are some more information about their output per-read methylation call files and columns for compatibility with patmat.py:  
nanopolish and f5c>=v0.7 produce the following columns and the CpG coordinates are zero-based and coordinates for both strands are based on positive strand (positions for the CpG from both strands are the same):
```
chromosome	strand	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_motifs	sequence
chr2	+	200000365	200000365	50152360-5abb-4e1f-9ce0-c08a49d65b57	3.91	-142.10	-146.01	1	1	GTGAACGCTTT
chr2	+	200000776	200000776	50152360-5abb-4e1f-9ce0-c08a49d65b57	-20.59	-243.72	-223.13	1	1	TAACTCGATTT
chr2	-	200000365	200000365	607a605c-f01b-4b02-a8d5-b4c8adb88e6b	4.93	-257.29	-262.22	1	1	GTGAACGCTTT
chr2	-	200000776	200000776	607a605c-f01b-4b02-a8d5-b4c8adb88e6b	-11.09	-225.59	-214.50	1	1	TAACTCGATTT
```
Megalodon per-read text methylation call output has the following columns and CpG coordinates are zero-based and coordinates of negative strand are 1 bp greater than positive strand (The methylation call file must be only for methylation. Do not use per-read methylation file that has multiple modification calls, e.g. 5mC and 5hmC):
```
read_id	chrm	strand	pos	mod_log_prob	can_log_prob	mod_base
56780a98-ccb3-41a5-8ed1-fc069412fc13    chr11   +       21488565        -0.9126647710800171     -0.5132502558405262     m
56780a98-ccb3-41a5-8ed1-fc069412fc13    chr11   +       21486004        -0.8042076826095581     -0.5931974211226271     m
2cc45d27-6084-49f1-b156-34501adc7651    chr11   -       21488566        -3.271272659301758      -0.03869726232984402    m
2cc45d27-6084-49f1-b156-34501adc7651    chr11   -       21486005        -4.3451995849609375     -0.013053750265459633   m
```
DeepSignal methylation call file has the following columns and CpG coordinates are zero-based and coordinates of negative strand are 1 bp greater than the positive strand:
```
chrom   pos     strand  pos_in_strand   readname        read_strand     prob_0  prob_1  called_label    k_mer
chr11	2669073	+	-1	19b5bd8e-0a50-449d-8dc1-ea2dc4e2fe2b	t	0.09740365	0.90259635	1	TACCCTGCCGTATCAGT
chr11	2669107	+	-1	19b5bd8e-0a50-449d-8dc1-ea2dc4e2fe2b	t	0.13432296	0.865677	1	ACTGGCTACGTGTGGCT
chr11	2669074	-	-1	12652f63-7676-4ad8-b7bf-af1aec4b282d	t	0.13398732	0.8660127	1	CACTGATACGGCAGGGT
chr11	2669108	-	-1	12652f63-7676-4ad8-b7bf-af1aec4b282d	t	0.12144542	0.87855464	1	GAGCCACACGTAGCCAG
```  

