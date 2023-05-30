#!/usr/bin/env Rscript

# Run like this: ./strandseq_phase.R -h

# Requires bcftools, R>=4.3.0, and the R packages devtools, BiocManager, InvertypeR, argparse, and BSgenome.Hsapiens.UCSC.hg38


#######################################################################################################
#######################################################################################################


# A modified error function so Rscript doesn't hide any errors that occur (https://renkun.me/2020/03/31/a-simple-way-to-show-stack-trace-on-error-in-r/)

err <- function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    sink(stderr())
    on.exit(sink(NULL))
    cat("Backtrace:\n")
    calls <- rev(calls[-length(calls)])
    for (i in seq_along(calls)) {
      cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
    }
  }
  if (!interactive()) {
    q(status = 1)
  }
}

options(error = err, warn = 1)

# A messy tryCatch routine that attempts to get the directory containing this script.
# This is then used as a default location for the blacklist and the inversion_list BED files, plus the soft_mask.
tryCatch(
    {
        initial_args <- commandArgs(trailingOnly = FALSE)
        script_dir <- dirname(sub("--file=", "", initial_args[grep("--file=", initial_args)]))
    },
    error=function(e) {
        print(e)
        print("The directory containing strandseq_phase.R could not be searched for --hard_mask, --soft_mask, and --inversion_list defaults")
    },
    warning=function(w) {
        print(w)
        print("The directory containing strandseq_phase.R could not be searched for --hard_mask, --soft_mask, and --inversion_list defaults")
    }
)

# The argparse library makes nice python-style command line arguments and flags
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description='Performs inversion-aware Strand-seq phasing of a VCF file of SNVs. Requires bcftools (samtools.github.io/bcftools/bcftools.html), 
    R>=4.3.0, and the R packages devtools (CRAN), BiocManager (CRAN), InvertypeR (GitHub: vincent-hanlon/InvertypeR), argparse (CRAN), and BSgenome.Hsapiens.UCSC.hg38 
    (Bioconductor).')

parser$add_argument("reference",
    type = "character",
    help = "Absolute path to the FASTA-format GRCh38 human reference genome.",
    metavar = "/path/to/GRCh38.fasta"
)

parser$add_argument("vcf",
    type = "character", 
    help = "Absolute path to a VCF file of SNVs to phase.",
    metavar = "/path/to/snps.vcf"
)

parser$add_argument("-p", "--paired",
    type = "logical", required = F, default = TRUE,
    help = "Are the Strand-seq reads paired end? Default: TRUE.",
    metavar = "TRUE or FALSE"
)

parser$add_argument("-i", "--input_folder",
    type = "character", required = F, default=".",
    help = "Absolute path to the directory containing good-quality Strand-seq libraries for your sample. Default: '.'.",
    metavar = "/path/to/input/"
)

parser$add_argument("-o", "--output_folder",
    type = "character", required = F, default=".",
    help = "Absolute path to a directory where output files should be written. Default: '.'.",
    metavar = "/path/to/output/"
)

parser$add_argument("-t", "--threads",
    type = "integer", required = F, default= 4,
    help = "The number of parallel threads to run on. Default: 4.",
    metavar = "integer"
)

parser$add_argument("-n", "--name",
    type = "character", required = F, default= "unknown",
    help = "The name of the sample. This will appear in the VCF file. E.g., 'HG005'. Default: 'unknown'.",
    metavar = "string"
)

parser$add_argument("--inversion_list",
    type = "character", required = F, default=file.path(script_dir, "hanlon_2021_BMCgenomics_augmented.bed"),
    help = "Absolute path to a BED file containing genomic intervals that might be inversions. This is typically a list from the literature
        so the file hanlon_2021_BMCgenomics_augmented.bed on the PatMat GitHub (originally from vincent-hanlon/InvertypeR) is a good start.
        Default: the file suggested above, if it is in the same directory as strandseq_phase.R.",
    metavar = "/path/to/BED"
)

parser$add_argument("--hard_mask",
    type = "character", required = F, default=file.path(script_dir, "hard_mask.GRCh38.humans.bed"),
    help = "Absolute path to a BED file containing regions with unreliable Strand-seq data. 
        The file hard_mask.GRCh38.humans.bed on the PatMat GitHub is a good start.
        Default: the file suggested above, if it is in the same directory as strandseq_phase.R.",
    metavar = "/path/to/BED"
)

parser$add_argument("--soft_mask",
    type = "character", required = F, default=file.path(script_dir, "soft_mask.bed"),
    help = "Absolute path to a BED file containing regions, like very large inversions, that
        occasionally interfere with composite file creation. Rarely really necessary (see InvertypeR documentation)
        Default: the file suggested above, if it is in the same directory as strandseq_phase.R, which contains
        the three largest autosomal inversions according to Porubsky et al. 2022, except for a very rare one on chr2.",
    metavar = "/path/to/BED"
)

parser$add_argument("--prior",
    type = "character", required = F, default = "0.9683,0.0158,0.0158", 
    help = "A comma-separated list of prior weights for inversion genotypes, without spaces. Only needs to be altered if --inversion_list is 
        not the default. E.g., 0.96,0.02,0.02. See InvertypeR for more details.",
    metavar = "comma-separated numbers"
)

chr <- "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

parser$add_argument("--chromosomes",
    type = "character", required = F, default = chr,
    help = "A comma-separated list of chromosome names, without spaces. Only SNVs on these chromosomes will be phased. 
        E.g., chr1,chr2,chr3. Default: the 22 autosomes.",
    metavar = "comma-separated strings"
)

args <- parser$parse_args()

stopifnot("Provide a valid VCF file with the --vcf flag" = file.exists(args$vcf))
stopifnot("Provide a valid FASTA reference genome with the --reference flag" = file.exists(args$reference))
stopifnot("Provide a valid path to a directory containing Strand-seq BAM files with the --input_folder flag" = file.exists(args$input_folder))
stopifnot("Provide a valid BED file of blacklisted regions with the --hard_mask flag" = file.exists(args$hard_mask))
stopifnot("Provide a valid BED file of putative inversions from the literature with the --inversion_list flag" = file.exists(args$inversion_list))

if(!file.exists(args$output_folder)){
    dir.create(args$output_folder)
}

suppressPackageStartupMessages(library("breakpointR"))
suppressPackageStartupMessages(library("StrandPhaseR"))
suppressPackageStartupMessages(library("invertyper"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))

# A function to reconcile (very crudely) inversions from two different sources: genotyping known coordinates and de novo discovery
combine_genotyped_discovered_inversions <- function(inversions) {

    if(nrow(inversions[[1]])==0){
        return(inversions[[2]])
    } else if(nrow(inversions[[2]])==0){
        return(inversions[[1]])
    } else {
        geno <- inversions[[1]][inversions[[1]]$end - inversions[[1]]$start + 1 > 10000 & inversions[[1]]$probability >= 0.95 &
            inversions[[1]]$genotype != 0 & inversions[[1]]$genotype != "0|0" & !inversions[[1]]$low_read_density, -10]

        disc <- inversions[[2]][inversions[[2]]$end - inversions[[2]]$start + 1 > 10000 & inversions[[2]]$probability >= 0.95 &
            inversions[[2]]$genotype != 0 & inversions[[2]]$genotype != "0|0" & !inversions[[2]]$low_read_density, -10]

        geno <- sort(GenomicRanges::makeGRangesFromDataFrame(geno, keep.extra.columns = T))
        disc <- sort(GenomicRanges::makeGRangesFromDataFrame(disc, keep.extra.columns = T))

        hits <- GenomicRanges::findOverlaps(geno, disc)
        overlaps <- GenomicRanges::pintersect(geno[S4Vectors::queryHits(hits)], disc[S4Vectors::subjectHits(hits)])
        percentOverlap1 <- GenomicRanges::width(overlaps) / GenomicRanges::width(disc[S4Vectors::subjectHits(hits)])
        percentOverlap2 <- GenomicRanges::width(overlaps) / GenomicRanges::width(geno[S4Vectors::queryHits(hits)])
        intersections <- hits[percentOverlap1 > 0.9 & percentOverlap2 > 0.9]
        new_geno <- geno[-S4Vectors::queryHits(intersections)]
        combined_inversions <- data.frame(sort(c(disc, new_geno)))

        return(combined_inversions)
    }
}

args$chromosomes <- unlist(strsplit(gsub(" ","",args$chromosomes),","))
args$prior <- as.numeric(unlist(strsplit(gsub(" ","",args$prior),",")))

message("\n       ##########     Inversions     ##########")

# calling inversions
# because this now uses an improved regions_to_genotype (better than the 2023 PofO paper), I have ...
# ... removed the smallest breakpointR windowsize = 40 and minReads = 15, which added ~ 6 hrs runtime on 12 CPU
inversions <- invertyper_pipeline(
    regions_to_genotype = args$inversion_list,
    prior = args$prior,
    adjust_method = "all",
    input_folder = args$input_folder,
    output_folder = file.path(args$output_folder, "invertyper"),
    haploid_chromosomes = NULL,
    vcf = args$vcf,
    paired_reads = args$paired,
    confidence = 0.95,
    hard_mask = args$hard_mask,
    soft_mask = args$soft_mask,
    chromosomes = args$chromosomes,
    numCPU = args$threads,
    save_composite_files = FALSE,
    write_browser_files = FALSE,
    discover_breakpointr_inversions = TRUE,
    breakpointr_prior = c(0.9, 0.05, 0.05),
    breakpointr_haploid_prior = c(0.9, 0.1),
    windowsize = c(120, 360),
    minReads = c(50, 50),
    background = 0.2,
    output_file = "inversions.txt"
)

# combining inversions from two sources and writing them to files
combined_inversions <- combine_genotyped_discovered_inversions(inversions)
write.table(combined_inversions, file = file.path(args$output_folder,"invertyper", paste0(args$name, ".inversions.above10kb.txt")), 
    sep = "\t", row.names = F, col.names = T, quote = F)
write.table(combined_inversions[, c(1:3)], file = file.path(args$output_folder, "invertyper", paste0(args$name, 
    ".inversions.above10kb.bed")), sep = "\t", row.names = F, col.names = F, quote = F)
message("\n       ##########     Phasing     ##########")
ptm <- startTimedMessage("\n       identifying WC regions for SNV phasing ...")

# finding Watson-Crick regions, which are suitable for phasing
invisible(pdf())
    suppressMessages(breakpointr(
        inputfolder = args$input_folder,
        outputfolder = file.path(args$output_folder, "BPR_output"),
        pairedEndReads = args$paired,
        numCPU = args$threads,
        windowsize = 8000000,
        binMethod = "size",
        chromosomes = args$chromosomes,
        background = 0.1,
        maskRegions = args$hard_mask
    ))
invisible(dev.off())

suppressMessages(exportRegions(file.path(args$output_folder, "BPR_output", "data"), file = "wc_regions.txt", collapseInversions = FALSE, 
    minRegionSize = 5000000, state = "wc"))

stopTimedMessage(ptm)
ptm <- startTimedMessage("\n       phasing SNVs (not inversion aware yet) ...")

# Phasing SNVs, (but getting inversions wrong)
suppressMessages(strandPhaseR(
    inputfolder = args$input_folder,
    outputfolder = file.path(args$output_folder, "SPR_output"),
    numCPU = args$threads,
    positions = args$vcf,
    WCregions = "wc_regions.txt",
    chromosomes = args$chromosomes,
    num.iterations = 3,
    exportVCF = args$name,
    bsGenome = "BSgenome.Hsapiens.UCSC.hg38",
    splitPhasedReads = TRUE,
    assume.biallelic = TRUE,
    pairedEndReads = args$paired
))

stopTimedMessage(ptm)
ptm <- startTimedMessage("       correcting phase at inversions ...")

# Correcting phasing within inversions
invisible(suppressMessages(correctInvertedRegionPhasing(
    outputfolder =  file.path(args$output_folder, "SPR_output", "VCFfiles"),
    inv.bed = file.path(args$output_folder, "invertyper", paste0(args$name, ".inversions.above10kb.bed")),
    strandphaseR.data = file.path(args$output_folder, "SPR_output", "data"),
    breakpointR.data = file.path(args$output_folder, "BPR_output", "data"),
    vcfs.files = file.path(args$output_folder, "SPR_output", "VCFfiles"),
    snv.positions = args$vcf,
    input.bams = args$input_folder,
    chromosomes = args$chromosomes,
    bsGenome = "BSgenome.Hsapiens.UCSC.hg38",
    ref.fasta = args$reference,
    recall.phased = TRUE,
    het.genotype = "lenient",
    pairedEndReads = args$paired,
    min.mapq = 10,
    background = 0.1,
    lookup.bp = 1000000,
    assume.biallelic = TRUE,
    lookup.blacklist = args$hard_mask
)))

stopTimedMessage(ptm)
ptm <- startTimedMessage("       concatenating VCF files ...")

# concatenating SPR's VCF files into a single file
# running bcftools from the command line seems easiest...
output_filename <- file.path(args$output_folder, paste0(args$name, ".strandseqphased.inv_aware.vcf"))
input_filenames <- file.path(args$output_folder, "SPR_output","VCFfiles","chr*INVcorr.vcf")
command <- paste0("bcftools concat ", input_filenames, " > ", output_filename)

invisible(system(command, ignore.stderr = TRUE))

stopTimedMessage(ptm)


