#! /usr/bin/env python3
# coding=utf-8

# Copyright (C) 2022  Vahid Akbari

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
PatMat: Parent-of-origin (Paternal and Maternal) resolved chromosome-scale haplotyping
"""

__author__ = "Vahid Akbari"
__email__ = "vakbari@bcgsc.ca"
__copyright__ = "Copyright (C) 2022, " + __author__
__license__ = "GPLv3"

import argparse
import os
import sys
import typing
import warnings

import pysam
import tabix

from patmat.core.allele_processing import process_variant
from patmat.core.haplotype import (
    add_reads_hap,
    build_read_dict_HP_temp,
    build_variant_dict_HP,
    update_per_var_info,
)
from patmat.core.methylation import PofO_dmr, out_freq_methbam, pofo_final_dict
from patmat.core.phase_processing import (
    alignment_writer,
    pofo_sv_write,
    write_sam_phase,
)
from patmat.core.subprocess_calls import bgzip_and_tabix, run_command
from patmat.core.variant_assignment import process_variant_assignments
from patmat.core.variant_processing import (
    process_strand_seq_vcf,
    write_per_read_variant_file,
)
from patmat.io.bam import get_chroms_from_bam
from patmat.io.file_utils import openfile
from patmat.io.vcf import get_chroms_from_vcf
from patmat.io.write import write_merged_dmr_regions, write_scores


def out_pofo_freq(hp_fre, chrom, output):
    """
    Outputs PofO assigned methylation frequency files.
    """
    if os.path.exists(hp_fre):
        with openfile(hp_fre) as file:
            next(file)
            for line in file:
                line = line.rstrip().split("\t")
                if line[0] == chrom:
                    output.write("\t".join(line) + "\n")
    else:
        with openfile(hp_fre + ".gz") as file:
            next(file)
            for line in file:
                line = line.rstrip().split("\t")
                if line[0] == chrom:
                    output.write("\t".join(line) + "\n")


def main(raw_arguments: typing.Optional[typing.List[str]] = None) -> None:
    args = _parse_arguments(raw_arguments or sys.argv[1:])
    """
    The main function that uses user's inputs and other functions to phase and
    assign PofO to variants and methylation.
    """
    hap_ratio = args.hapratio
    min_variant = args.min_variant
    min_cg = args.min_cg
    bam_file = os.path.abspath(args.bam)
    processes = args.processes
    dss_processes = args.dss_processes or processes
    cpg_difference = args.cpg_difference
    chunk = args.chunk_size
    min_read_reassignment = args.min_read_number
    out_prefix = os.path.abspath(args.output)
    reference = os.path.abspath(args.reference)
    vcf = os.path.abspath(args.vcf)
    known_dmr = os.path.abspath(args.known_dmr)

    if not os.path.isfile(vcf + ".tbi"):
        raise Exception("It seems that vcf file " "is not index by tabix.")
    vcf_tb = tabix.open(vcf)

    if args.black_list is not None:
        sites_to_ignore = build_sites_to_ignore(vcf, args.black_list)

    chroms = get_chroms_from_vcf(vcf)
    bam_choms = get_chroms_from_bam(bam_file)

    reads_hap = dict()
    re_assignment_vars = dict()
    for chrom in sorted(chroms.keys()):
        print("#############  Processing chromosome {}  #############".format(chrom))

        final_dict, strand_phased_vars, phase_block_stat, blocks_dict = (
            process_strand_seq_vcf(
                vcf,
                chrom,
                args.strand_vcf,
                args.phased,
                args.include_all_variants,
                args.ignore_blocks_single,
            )
        )

        if strand_phased_vars == 0:
            warnings.warn(f"No phased strand-seq variant for {chrom}, skipping it.")
            continue

        read_info_file = out_prefix + "_temp_VarReadinfo_" + chrom + ".tsv"
        write_per_read_variant_file(
            final_dict,
            bam_file,
            chunk,
            processes,
            read_info_file,
            args.mapping_quality,
            args.include_supplementary,
        )
        read_dict_HP_temp = build_read_dict_HP_temp(read_info_file)
        if not read_dict_HP_temp:
            warnings.warn("No phased read for {}." " Skipping it.".format(chrom))
            continue

        per_var_info, read_dict_HP_temp_reass = update_per_var_info(
            hap_ratio,
            min_variant,
            read_info_file,
            read_dict_HP_temp,
        )

        add_reads_hap(hap_ratio, min_variant, reads_hap, read_dict_HP_temp_reass)

        variant_dict_HP = build_variant_dict_HP(reads_hap, read_info_file, per_var_info)

        records_chrom = vcf_tb.query(chrom, 0, chroms[chrom] + 1)
        for line in records_chrom:
            format_values = dict(zip(line[8].split(":"), line[9].split(":")))
            if "PS" in format_values:
                del format_values["PS"]

            if not args.include_all_variants and line[6] not in ["PASS", "."]:
                continue
            if format_values["GT"].startswith(
                ("0/1", "1/0", "0|1", "1|0", "1/2", "1|2", "2/1", "2|1")
            ):
                re_assignment_vars[tuple(line[0:2])] = process_variant(
                    line,
                    format_values,
                    variant_dict_HP,
                    per_var_info,
                    min_read_reassignment,
                    phase_block_stat,
                    blocks_dict,
                )
        os.remove(read_info_file)

        per_var_info.clear()

    temp_DMR_file = out_prefix + "_temp_knownDMR.tsv"
    write_merged_dmr_regions(known_dmr, temp_DMR_file)

    temp_at_DMRs_file = out_prefix + "_TempAtDMRs.sam.gz"
    run_command(
        "samtools view -h --remove-tag HP,PS -@ {} "
        "--regions-file {} {}"
        " | sed '/^@PG/d' | gzip > {}".format(
            processes, temp_DMR_file, bam_file, temp_at_DMRs_file
        )
    )
    run_command("rm {}*".format(out_prefix + "_temp"))

    temp_non_pofo_dmr_file = out_prefix + "_Temp-NonPofO_dmr.sam.gz"
    write_sam_phase(
        temp_at_DMRs_file,
        temp_non_pofo_dmr_file,
        reads_hap,
        args.mapping_quality,
        args.include_supplementary,
    )

    temp_non_fomo_dmr_bam_file = out_prefix + "_Temp-NonPofO_dmr.bam"
    run_command(
        "gunzip -c {0} | samtools sort -@ {1} -o {2} && "
        "samtools index -@ {1} {2}".format(
            temp_non_pofo_dmr_file, processes, temp_non_fomo_dmr_bam_file
        ),
    )
    out_freqhp1 = out_prefix + "_Temp_NonPofO_HP1-HP2_MethylationHP1.tsv"
    out_freqhp2 = out_prefix + "_Temp_NonPofO_HP1-HP2_MethylationHP2.tsv"
    out_freq_methbam(
        out_prefix,
        processes,
        reference,
        args.pacbio,
    )

    run_command(
        "{} {} {} {} {} {} {} {} {} {} {}".format(
            "Rscript",
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "DMA_UsingDSS.R"),
            out_freqhp1,
            out_freqhp2,
            out_prefix + "_Temp",
            args.equal_disp,
            args.smoothing_flag,
            args.smoothing_span,
            args.delta_cutoff,
            args.pvalue,
            dss_processes,
        ),
    )
    bgzip_and_tabix(out_prefix + "_Temp_DMLtest.tsv")
    bgzip_and_tabix(out_prefix + "_Temp_callDML.tsv")

    chrom_hp_origin_count = PofO_dmr(known_dmr, out_prefix, min_cg, cpg_difference)
    chrom_hp_origin = pofo_final_dict(chrom_hp_origin_count, args.min_pofo_score)

    print("################## Assigning PofO to Variants ##################")
    info_out_dict = process_variant_assignments(
        vcf, out_prefix, re_assignment_vars, chrom_hp_origin
    )

    if args.sv_vcf is not None:
        print("################## Assigning PofO to SVs ##################")
        all_sv_files = args.sv_vcf
        for sv_file in all_sv_files:
            pofo_sv_write(
                sv_file,
                out_prefix,
                chrom_hp_origin,
                reads_hap,
                min_read_reassignment,
                args.include_all_variants,
                hap_ratio,
            )

    print("############ Preparing PofO Tagged Alignment File #############")
    pofo_tagged_cram_file = out_prefix + "_PofO_Tagged.cram"
    with pysam.AlignmentFile(bam_file, "rb") as infile_bam:
        cramHeader = infile_bam.header.to_dict()
        if "PG" in cramHeader:
            cramHeader = cramHeader.pop("PG")
        with pysam.AlignmentFile(
            pofo_tagged_cram_file,
            "wc",
            template=infile_bam,
            header=cramHeader,
            reference_filename=reference,
        ) as outfile_bam:
            for chrom in list(bam_choms) + ["*"]:
                alignment_writer(
                    infile_bam,
                    chrom,
                    reads_hap,
                    chrom_hp_origin,
                    outfile_bam,
                    args.mapping_quality,
                    args.include_supplementary,
                )

    run_command(
        "samtools sort -@ {0} {1} -o {1} && samtools index -@ {0} -c {1}"
        "".format(processes, pofo_tagged_cram_file),
    )
    run_command("rm {}*".format(out_prefix + "_Temp"))

    write_scores(out_prefix, chrom_hp_origin)

    if not args.include_all_variants:
        print(
            'Per chromosome info for the variants with "PASS"'
            ' or "." in FILTER column in the input vcf file:'
        )
    else:
        print("Per chromosome info for the variants in the input vcf file:")
    print(
        "chrom\tall_het_snvs\tall_het_indels\t"
        "pofo_assigned_het_snvs\tpofo_assigned_het_indels"
    )
    for key, val in info_out_dict.items():
        print(
            "\t".join(
                map(
                    str,
                    [
                        key,
                        val["all_het_snvs"],
                        val["all_het_indels"],
                        val["pofo_het_snvs"],
                        val["pofo_het_indels"],
                    ],
                )
            )
        )


def build_sites_to_ignore(vcf, black_list):
    sites_to_ignore = set()
    if not os.path.isfile(vcf + ".tbi"):
        raise Exception(
            "Black list is given but it seems that"
            " the vcf file is not indexed (file ends with"
            " .tbi was not found)."
        )
    tb_vcf = tabix.open(vcf)
    black_list_file = os.path.abspath(black_list)
    with openfile(black_list_file) as bl:
        for line in bl:
            line = line.rstrip().split("\t")
            try:
                records = tb_vcf.query(line[0], int(line[1]), int(line[2]) + 1)
            except:
                warnings.warn(
                    "{}:{}-{} region from black list does not exist in the "
                    "vcf file. Skipping it.".format(line[0], line[1], line[2])
                )
                records = "NA"
            if records != "NA":
                for record in records:
                    sites_to_ignore.add((record[0], str(int(record[1]) - 1)))
    return sites_to_ignore


def _parse_arguments(raw_arguments: typing.List[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        raw_arguments (List[str]): list of arguments

    Returns:
        argparse.Namespace: namespace object with parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog="patmat",
        add_help=False,
        description="Phasing reads and Methylation "
        "using strand-seq and nanopore to determine "
        "PofO of each homologous chromosome "
        "in a single sample.",
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "--bam",
        "-b",
        action="store",
        type=str,
        required=True,
        help="The path to the coordinate sorted bam file with methylation" " tag.",
    )
    required.add_argument(
        "--output",
        "-o",
        action="store",
        type=str,
        required=True,
        help=(
            "The path to directory and prefix to save "
            "files. e.g path/to/directory/prefix"
        ),
    )
    required.add_argument(
        "--vcf",
        "-v",
        action="store",
        type=str,
        required=True,
        default=None,
        help="The path to the vcf file. If the input vcf is phased (e.g. using whatshap or longphase) "
        "you can select the --phased option to also use phase blocks. See --phased for more details.",
    )
    required.add_argument(
        "--strand_vcf",
        "-stv",
        action="store",
        type=str,
        required=True,
        help="The path to the chromosome-scale phased vcf file."
        " This is the input vcf file that has been phased "
        "using strand-seq data.",
    )
    required.add_argument(
        "--reference",
        "-ref",
        action="store",
        type=str,
        required=False,
        help=(
            "If you have given a bam file with methylation tag"
            " for the --tool_and_callthresh option, then you "
            "must also give the path to the reference file. "
            "File must be indexed using samtools faidx."
        ),
    )
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument(
        "--min_pofo_score",
        "-mps",
        action="store",
        type=float,
        required=False,
        default=0.55,
        help=(
            "0-1. Threshold for chromosome PofO score. This is "
            "the minimum PofO score of a chromosome to assign PofO "
            "to a haplotype. See the github page on how PofO score "
            "of a chromosome is calculated. Default is 0.55"
        ),
    )
    optional.add_argument(
        "--pacbio",
        "-pb",
        action="store_true",
        required=False,
        help="Select this if the reads are from PacBio HiFi.",
    )

    optional.add_argument(
        "--known_dmr",
        "-kd",
        action="store",
        type=str,
        required=False,
        default=os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "Imprinted_DMR_List_V1.GRCh38.tsv",
        ),
        help="The path to the input file for known imprinted DMRs. "
        "File must have the following information in the "
        "following column order: "
        "chromosome\tstart\tend\tMethylatedAlleleOrigin "
        "where the methylated allele origin must be either "
        "maternal or paternal (First row must be header). "
        "By default, we use hg38 reference iDMR list in repo's patmat directory."
        " We have also provided the iDMR list for t2t (hs1) reference"
        " in the  repo's patmat directory, in case you are using t2t "
        "reference provide the path to that file.",
    )
    optional.add_argument(
        "--phased",
        "-ph",
        action="store_true",
        required=False,
        help=(
            "Select this option if your input vcf is a phased vcf file"
            " using long-reads using whatshap or longphase. "
            "In this case phased-blocks and phased variants "
            "inside them will be used and strand-seq phased data"
            " are used to correct the switches accross phase blocks. "
            "This can be useful when the chromosome-scale"
            " phased variants from strand-seq are very sparse. "
            "If selected, vcf file must be indexed using tabix."
            " Note that phased blocks will be extract from vcf file and"
            ' assumption is that the number after last ":" sign in '
            "column 10 is the block ID, as with vcfs phased by "
            "whatshap or longphase."
        ),
    )
    optional.add_argument(
        "--ignore_blocks_single",
        "-ibs",
        action="store_true",
        required=False,
        help=(
            "When --phased is selected, strand-seq phased variants "
            "at phase blocks will be used for switch correction "
            "by a two variant overlap window. However, some phase "
            "blocks have only a single variant from strand-seq."
            " By default those blocks will also be considered. "
            "Select this option if you wish to ignore them."
        ),
    )

    optional.add_argument(
        "--sv_vcf",
        "-sv_vcf",
        action="store",
        type=str,
        nargs="+",
        required=False,
        default=None,
        help=(
            "Path to the Structural variation (SV) call"
            " file if you wish to also add PofO to SVs."
            " if multiple files are gived, they must be separated by"
            " space and give the absolute path to each file."
            " File(s) must include read names (e.g.RNAMES=read1,"
            "read2,read3) in the 8th column of vcf ."
        ),
    )
    optional.add_argument(
        "--black_list",
        "-bl",
        action="store",
        type=str,
        required=False,
        default=None,
        help="List of regions to ignore phased variants at them."
        " Three first columns must be chromosome\tstart\tend."
        " If a black list is given the vcf file must be indexed using tabix.",
    )
    optional.add_argument(
        "--hapratio",
        "-hr",
        action="store",
        type=float,
        required=False,
        default=0.75,
        help=(
            "0-1. For phasing reads, minimum ratio of phased variants "
            "a read must have from a haplotype to assign it to that haplotype. "
            "For phasing SVs, minimum ratio of phased reads on a haplotype. "
            "Default is 0.75."
        ),
    )
    optional.add_argument(
        "--mapping_quality",
        "-mq",
        action="store",
        type=int,
        required=False,
        default=20,
        help=(
            "An integer value to specify threshold for "
            "filtering reads based on mapping quality for "
            "PofO assignment to variants. "
            "Default is >=10"
        ),
    )
    optional.add_argument(
        "--min_variant",
        "-mv",
        action="store",
        type=int,
        required=False,
        default=1,
        help=(
            "Minimum number of phased variants a read must have "
            "to be considered during variant rephasing."
            ". Default= 1."
        ),
    )
    optional.add_argument(
        "--min_read_number",
        "-mr",
        action="store",
        type=int,
        required=False,
        default=2,
        help=(
            "Minimum number of reads to support a variant"
            " to assign to each haplotype. Default= 2"
        ),
    )
    optional.add_argument(
        "--min_cg",
        "-mcg",
        action="store",
        type=int,
        required=False,
        default=5,
        help=(
            "Minimum number of CpGs an iDMR must have to "
            " consider it for PofO assignment. Default is 5."
        ),
    )
    optional.add_argument(
        "--cpg_difference",
        "-cd",
        action="store",
        type=float,
        required=False,
        default=0.1,
        help=(
            "Minimum cut off for the fraction of CpGs between haplotypes "
            "must be differentially methylated at an iDMR to "
            "consider it for PofO assignment. Default is 0.1."
        ),
    )
    optional.add_argument(
        "--include_all_variants",
        "-iav",
        action="store_true",
        required=False,
        help='By default, only variants that have "PASS" or "." '
        " in the FILTER column of the input vcf file will be used"
        " during phasing and PofO assignment. Select this flag "
        "if you want to use all the variants.",
    )
    optional.add_argument(
        "--include_supplementary",
        "-is",
        action="store_true",
        required=False,
        help="Include supplementary reads.",
    )
    optional.add_argument(
        "--processes",
        "-p",
        action="store",
        type=int,
        required=False,
        default=4,
        help="Number of parallel processes. Default is 4.",
    )
    optional.add_argument(
        "--chunk_size",
        "-cs",
        action="store",
        type=int,
        required=False,
        default=500,
        help=("Chunk per process. Default is 500"),
    )
    optional = parser.add_argument_group(
        "Optional arguments. The following options "
        "are DSS options for differential methylation"
        " analysis to find differentially methylated "
        "CpGs between haplotypes"
    )
    optional.add_argument(
        "--delta_cutoff",
        "-dc",
        action="store",
        type=float,
        default=0.075,
        required=False,
        help=(
            "0-1. A threshold for defining differentially "
            "methylated loci (DML) or CpGs. In DML testing"
            " procedure, hypothesis test that the two groups "
            "means are equal is conducted at each CpG site. "
            "Here if delta is specified, the function will "
            "compute the posterior probability that the "
            "difference of the means are greater than delta,"
            " and then call DML based on that. Default is 0.075."
        ),
    )
    optional.add_argument(
        "--pvalue",
        "-pv",
        action="store",
        type=float,
        required=False,
        default=0.001,
        help=(
            "0-1. When delta is not specified, this is the "
            "threshold of p-value for defining DML and "
            "loci with p-value less than this threshold "
            "will be deemed DMLs. When delta is specified, "
            "CpG sites with posterior probability greater than"
            " 1-pvalue_threshold are deemed DML. Default is 0.001"
        ),
    )
    optional.add_argument(
        "--smoothing_span",
        "-sms",
        action="store",
        type=int,
        default=500,
        required=False,
        help=("The size of smoothing window, in " "basepairs. Default is 500."),
    )
    optional.add_argument(
        "--smoothing_flag",
        "-smf",
        action="store",
        type=str,
        default="TRUE",
        required=False,
        help=(
            "TRUE/FALSE. A flag to indicate whether to apply "
            "smoothing in estimating mean methylation levels."
            " For more instructions see the DSS R package guide. "
            "Default is TRUE."
        ),
    )
    optional.add_argument(
        "--equal_disp",
        "-ed",
        action="store",
        type=str,
        default="FALSE",
        required=False,
        help=(
            "TRUE/FALSE. A flag to indicate whether the "
            "dispersion in two groups are deemed equal or not. "
            "For more instructions see the DSS R package guide. "
            "Default is FALSE Because there is no biological"
            " replicate here, you should specify either "
            "equal_disp TRUE or smoothing_flag TRUE. "
            "Do not specify both as FALSE."
        ),
    )
    optional.add_argument(
        "--dss_processes",
        "-dp",
        action="store",
        type=int,
        required=False,
        help="Number of parallel processes use for DSS "
        "differential methylation analysis. If not given, it will "
        "be the same as --processes option. Differential methylation "
        " analysis usually requires high memory and if there"
        " is not enough memory, specify less number processes"
        " using this flag for DSS to allocate available memory "
        "for less processes.",
    )
    optional = parser.add_argument_group("Help and version options")
    optional.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.4.0",
        help="Print program's version and exit",
    )
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Print this help and exit.",
    )
    parsed_arguments = parser.parse_args(raw_arguments)
    return parsed_arguments


if __name__ == "__main__":
    main()
