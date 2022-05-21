![](Haplotype.png)
# PatMat:  
This workflow enables simultanous chromosome-scale haplotyping and parent-of-origin detection using a combination of nanopore sequencing and Strand-seq.
We will used nanopore-detected variants and their long-range phasing from Strand-seq to detect chromosome-scale haplotypes. We then use DNA methylation information at known imprinted regions to detect parent-of-origin.  

To run this workflow you will need the following third-party tools:  
**Guppy**: For basecalling nanopore reads.  
**[Minimap2](https://github.com/lh3/minimap2)**: To align nanopore reads to the reference genome.  
**[Nanopolish](https://github.com/jts/nanopolish)**: To call DNA methylation from nanopore data.  
**[NanoMethPhase](https://github.com/vahidAK/NanoMethPhase)**: To process methylation call results from nanopolish.  
**[Clair3](https://github.com/HKU-BAL/Clair3)**: To call variants from aligned nanopore reads.  

The workflow was developed using the above tools. However, you may use alternatives to each tool.  

Finally you need to use our tool in this repository "**PatMat.py**" to detect chromosome-scale parent-of-origin resolved haplotypes.  
Table of Contents
=================

* **[Full Tutorial](https://github.com/vahidAK/PatMat/blob/master/README.md#full-tutorial)**
  * [Basecalling from nanopore data](https://github.com/vahidAK/PatMat/blob/master/README.md#1--Basecalling-from-nanopore-data)
  * [Mapping nanopore basecalled reads](https://github.com/vahidAK/PatMat/blob/master/README.md#1--Mapping-nanopore-basecalled-reads)
  * [Methylation Calling from nanopore data](https://github.com/vahidAK/PatMat/blob/master/README.md#1--Methylation-Calling-from-nanopore-data)
  * [Variant Calling from nanopore data](https://github.com/vahidAK/PatMat/blob/master/README.md#2--Variant-Calling-from-nanopore-data)
  * [Phasing variants using Strand-sq data](https://github.com/vahidAK/PatMat/blob/master/README.md#3--Phasing-variants-using-Strand\-sq-data)
  * [Parent-of-origin detection](https://github.com/vahidAK/PatMat/blob/master/README.md#4--Parent\-of\-origin-detection)
  
# Full Tutorial


## 1- Basecalling from nanopore data
We use Oxford Nanopore Technologies' basecaller, guppy, for translating raw signals to DNA sequence.
```
guppy_basecaller --input_path <Path to fast5 directory> \
  --save_path <path to output directory> \
  --config <appropriate configuration file (e.g. dna_r9.4.1_450bps_sup_prom.cfg)> \
  --device <GPU devices to be used (e,g. cuda:0 cuda:1)> \
  --trim_strategy dna
```
After basecalling you need to merge all the fastq files to a single file.

## 2- Mapping nanopore basecalled reads
Here we use [minimap2](https://github.com/lh3/minimap2) to align nanopore reads to reference genome. You may use other tools such as [Winnowmap](https://github.com/marbl/Winnowmap).
```
minimap2 -ax map-ont --MD -L -t <# of threads> \
  /path/to/reference.fa \
  /path/to/reads.fastq > /path/to/aligned_reads.sam 
```
After alignment was complete you need to sort and index alignment file.
```
samtools sort -@ <# of threads> /path/to/aligned_reads.sam  -o /path/to/aligned_reads.bam 
samtools index -@ <# of threads> /path/to/aligned_reads.bam
```

## 3- Methylation Calling from nanopore data  
Here we use [nanopolish](https://github.com/jts/nanopolish) for methylation calling however you may use [megalodon](https://github.com/nanoporetech/megalodon) or [DeepSignal](https://github.com/bioinfomaticsCSU/deepsignal). 

### 3-1 indexing fastq file and fast5 files:

NOTE: Fastqs must be merged to a single file

```
nanopolish index -d /path/to/reads.fastq
```

### 3-2 Methylation calling for CpG from each read:

```
nanopolish call-methylation \
  -t <number_of_threads> -q cpg \
  -r /path/to/reads.fastq \
  -b /path/to/aligned_reads.bam \
  -g /path/to/reference.fa > /path/to/MethylationCall.tsv
```

For the full tutorial please refer to
[Nanopolish](https://github.com/jts/nanopolish) page on GitHub.
### 3-3 Pre-processing methylation call file
We then need to pre-process methylation call file from nanopolish using [NanoMethPhase](https://github.com/vahidAK/NanoMethPhase) methyl_call_processor module.
```
nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > MethylationCall.bed.gz && tabix -p bed MethylationCall.bed.gz
```
## 4- Variant Calling from nanopore data

Here use [Clair3](https://github.com/HKU-BAL/Clair3) to call variants. However, you may call variants with other
tools such as [deepvariant](https://github.com/google/deepvariant).

```
run_clair3.sh --bam_fn=/path/to/aligned_reads.bam \
  --ref_fn=/path/to/reference.fa \
  --output=/path/to/output/directory \
  --threads=<# of threads> --platform=ont \
  --model_path=/path/to/model/ont_guppy5_r941_sup_g5014
```
After variant calling the results will be in merge_output.vcf.gz file in the output directory. You then need to extract high uality variants:  
```
gunzip -c /path/to/output/directory/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /path/to/output/Passed_Clair3_Variants.vcf
```  

## 5- Phasing variants using Strand-sq data
  

## 5- Parent-of-origin detection


