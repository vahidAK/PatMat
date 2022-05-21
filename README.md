![](Haplotype.png)
# PatMat:
  
Table of Contents
=================

* **[Full Tutorial](https://github.com/vahidAK/PatMat/blob/master/README.md#full-tutorial)**
  * [Methylation Calling from nanopore data](https://github.com/vahidAK/PatMat/blob/master/README.md#1--Methylation-Calling-from-nanopore-data)
  * [Variant Calling from nanopore data](https://github.com/vahidAK/PatMat/blob/master/README.md#2--Variant-Calling-from-nanopore-data)
  * [Phasing variants using Strand-sq data](https://github.com/vahidAK/PatMat/blob/master/README.md#3--Phasing-variants-using-Strand\-sq-data)
  * [Parent-of-origin detection](https://github.com/vahidAK/PatMat/blob/master/README.md#4--Parent\-of\-origin-detection)
  
# Full Tutorial

In order to get the phased methylome you also need the following third-party
software:

[Nanopolish](https://github.com/jts/nanopolish) : To call CpG methylation.

[Clair](https://github.com/HKU-BAL/Clair) or other variant callers: To call
variants for your sample. Alternatively, you might already have variant calling
data for example from Illumina sequencing.

[WhatsHap](https://github.com/whatshap/whatshap): To phase single nucleotide
variants.

## 1- Methylation Calling from nanopore data

### 1-1 indexing fastq file and fast5 files:

NOTE: Fastqs must be merged to a single file

```
nanopolish index -d /path/to/fast5s_directory/.fastq
```

### 1-2 Methylation calling for CpG from each read:

```
nanopolish call-methylation -t <number_of_threads> -q cpg -r /path/to/fastq_fromstep-1/fastq.fastq -b /path/to/sorted_and_indexed/bam.bam -g /path/to/reference.fa > /path/to/MethylationCall.tsv
```

For the full tutorial please refer to
[Nanopolish](https://github.com/jts/nanopolish) page on GitHub.

## 2- Variant Calling from nanopore data

We have used Clair to call variants. However, you may call variants with other
tools or your variant data may come from Illumina or other methods.

You can call variants for each chromosome using the following command and then
concatenate all files:

```
for i in chr{1..22} chrX chrY; do callVarBam --chkpnt_fn <path_to_model_file> --ref_fn <reference_genome.fa> --bam_fn <sorted_indexed.bam> --ctgName $i --sampleName <your_sample_name> --call_fn $i".vcf" --threshold 0.2 --samtools <path_to_executable_samtools_software> --pypy <path_to_executable_pypy > --threads <number_of_threads>
```

For the full tutorial please refer to [Clair](https://github.com/HKU-BAL/Clair)
page on GitHub.

After variant calling, you can select only SNVs which will be used for phasing:
```
awk '$4 != "." && $5 != "." && length($4) == 1 && length($5) == 1 && $6 > <the_variant_calling_quality_threshold>' variants.vcf > HighQualitySNVs.vcf
```

If you are calling variants from low coverage nanopore data (<30x) using Clair, you can also use our other tool [SNVoter](https://github.com/vahidAK/SNVoter) to improve SNV detection.

## 3- Phasing variants using Strand-sq data

  
**NOTE:** Currently, NanoMethPhase requires a single sample vcf file in which phase information of SNVs in 10th column indicated by "|" (e.g. 0|1 or 1|0). For more information about the input vcf file please read issue [#1](https://github.com/vahidAK/NanoMethPhase/issues/1).  
## 4- Parent-of-origin detection

### 4-1 First you need to process methylation call file using [NanoMethPhase](https://github.com/vahidAK/NanoMethPhase) methyl_call_processor module.

```
nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > MethylationCall.bed.gz && tabix -p bed MethylationCall.bed.gz
```
