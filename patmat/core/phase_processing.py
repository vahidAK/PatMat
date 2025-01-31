import gzip
import os
import warnings
from collections import defaultdict
from typing import DefaultDict, Dict, List, Optional, Set, TextIO, Tuple, Union

import pysam

from patmat.io.file_utils import openfile


def get_block(
    vcf: str, chrom: str
) -> Tuple[List[Tuple[str, int, int]], Dict[Tuple[str, str], Tuple[str, str]]]:
    """Extract phased blocks from VCF file for a given chromosome.

    Args:
        vcf: Path to VCF file
        chrom: Chromosome to process

    Returns:
        Tuple containing:
            - List of (chrom, start, end) tuples for each block
            - Dictionary mapping (chrom, block_id) to (start, end) strings
    """

    def extract_blocks(line: List[str]) -> Tuple[Tuple[str, str], int]:
        """Extract block information from a VCF line."""
        block_id = line[9].split(":")[-1]
        position = int(line[1])
        return ((line[0], block_id), position)

    blocks_dict: DefaultDict[Tuple[str, str], Set[int]] = defaultdict(set)

    # Extract block positions
    with openfile(vcf) as vp:
        for line in vp:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            if line[0] != chrom or "|" not in line[9].split(":")[0]:
                continue
            key, pos = extract_blocks(line)
            blocks_dict[key].add(pos)

    # Build final block lists
    final: List[Tuple[str, int, int]] = []
    formatted_blocks: Dict[Tuple[str, str], Tuple[str, str]] = {}

    for key, positions in blocks_dict.items():
        sorted_pos = sorted(positions)
        final.append((key[0], sorted_pos[0], sorted_pos[-1]))
        formatted_blocks[key] = (str(sorted_pos[0]), str(sorted_pos[-1]))

    return final, formatted_blocks


def get_sv_pofo(
    feed_list: List[Tuple[str, int, Set[str]]], alignment_file: str
) -> Tuple[DefaultDict[Tuple[str, int], int], DefaultDict[Tuple[str, int], int]]:
    """Map reads to haplotypes for structural variants.

    Args:
        feed_list: List of (chrom, position, read_ids) tuples
        alignment_file: Path to BAM/CRAM file

    Returns:
        Tuple of two defaultdicts mapping (chrom, pos) to counts for each haplotype
    """

    def process_pileup_read(
        pileupread: pysam.PileupRead,
        sv_reads: Set[str],
        sv_hp_dict1: DefaultDict[Tuple[str, int], int],
        sv_hp_dict2: DefaultDict[Tuple[str, int], int],
        position: Tuple[str, int],
    ) -> None:
        """Process a single pileup read and update haplotype counts."""
        read_id = pileupread.alignment.query_name
        try:
            hp_tag = pileupread.alignment.get_tag("HP")
        except:
            return

        if read_id in sv_reads:
            if hp_tag == 1:
                sv_hp_dict1[position] += 1
                sv_hp_dict2[position] += 0
            elif hp_tag == 2:
                sv_hp_dict1[position] += 0
                sv_hp_dict2[position] += 1

    sv_hp_dict1: DefaultDict[Tuple[str, int], int] = defaultdict(int)
    sv_hp_dict2: DefaultDict[Tuple[str, int], int] = defaultdict(int)

    samfile = pysam.AlignmentFile(alignment_file, "rb")

    for chrom, position, sv_reads in feed_list:
        try:
            sam_pileup = samfile.pileup(chrom, position, position + 1, truncate=True)
        except:
            warnings.warn(
                f"Variant {chrom} {position + 1} not found or no mapped "
                f"reads in alignment file. Check if correct BAM file is "
                f"given or BAM is indexed and not corrupted. Skipping it."
            )
            continue

        for pileupcolumn in sam_pileup:
            pileupcolumn.set_min_base_quality(0)
            if pileupcolumn.pos == position:
                for pileupread in pileupcolumn.pileups:
                    process_pileup_read(
                        pileupread,
                        sv_reads,
                        sv_hp_dict1,
                        sv_hp_dict2,
                        (chrom, position),
                    )

    return sv_hp_dict1, sv_hp_dict2


def count_reads_haplotypes(
    chromosome: str, reads_hap: Dict[Tuple[str, str], int], read_column: str
) -> DefaultDict[int, int]:
    """Count reads for each haplotype from the RNAMES column.

    Args:
        chromosome: Chromosome name
        reads_hap: Dictionary mapping (chrom, read_id) to haplotype assignment (1 or 2)
        read_column: VCF column containing read names (RNAMES=read1,read2,...)

    Returns:
        DefaultDict mapping haplotype (1 or 2) to read count
    """
    hp_counts = defaultdict(int)
    for read_ID in read_column.split("RNAMES=")[1].split(";")[0].split(","):
        hap = reads_hap.get((chromosome, read_ID))
        if hap:
            hp_counts[hap] += 1
    return hp_counts


def determine_haplotype_from_counts(
    hp_counts: DefaultDict[int, int], hapratio: float, min_read_reassignment: int
) -> Optional[str]:
    """Determine dominant haplotype based on read counts and thresholds.

    Args:
        hp_counts: Dictionary mapping haplotype (1 or 2) to read count
        hapratio: Minimum ratio required for haplotype assignment
        min_read_reassignment: Minimum number of reads required for reassignment

    Returns:
        "HP1", "HP2", or None if no haplotype is dominant
    """
    for hap in [1, 2]:
        if (
            hp_counts[hap] >= min_read_reassignment
            and hp_counts[hap] / (hp_counts[1] + hp_counts[2]) >= hapratio
        ):
            return f"HP{hap}"
    return None


def write_sv_line(
    line: List[str],
    format_values: Dict[str, str],
    additional_values: List[str],
    replace_patterns: List[Tuple[str, str]],
    sv_assignment_file: TextIO,
    sv_assignment_info_file: TextIO,
) -> None:
    """Write a single structural variant line to output files.

    Args:
        line: List of VCF fields
        format_values: Dictionary of FORMAT field values
        additional_values: List of additional values to append
        replace_patterns: List of string replacement tuples (old, new)
        sv_assignment_file: Main VCF output file handle
        sv_assignment_info_file: Info file handle for extra statistics
    """
    value_string = ":".join(format_values.values())
    for pattern in replace_patterns:
        value_string = value_string.replace(pattern[0], pattern[1])

    line_out = (
        line[0:8]
        + [":".join(format_values.keys())]
        + [value_string]
        + additional_values
    )
    sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
    sv_assignment_info_file.write("\t".join(line_out) + "\n")


def pofo_sv_write(
    sv_file: str,
    out: str,
    chrom_hp_origin: Dict[str, Dict[str, List[Union[str, int, float]]]],
    reads_hap: Dict[Tuple[str, str], int],
    min_read_reassignment: int,
    include_all_variants: bool,
    hapratio: float,
) -> None:
    """Write parent-of-origin assignments for structural variants.

    Process structural variants from VCF file and write parent-of-origin assignments
    based on haplotype read counts and chromosome assignments.

    Args:
        sv_file: Input structural variants VCF file path
        out: Output file prefix path
        chrom_hp_origin: Dictionary mapping chromosomes to haplotype parent-of-origin assignments
                        and statistics. Structure: {chrom: {"HP1": [origin, stats...], "HP2": [origin, stats...]}}
        reads_hap: Dictionary mapping (chromosome, read_id) tuples to haplotype assignment (1 or 2)
        min_read_reassignment: Minimum number of reads required to reassign a variant
        include_all_variants: Whether to process variants that don't pass filters
        hapratio: Minimum ratio of reads required to assign a haplotype
    """

    # Start
    sv_assignment_file = open(
        out + "_" + os.path.basename(sv_file) + "_PofO_Assignment_SVs.vcf", "w"
    )
    sv_assignment_info_file = open(
        out + "_" + os.path.basename(sv_file) + "_Variant_Assignment_SVs_info.tsv", "w"
    )
    with openfile(sv_file) as vf:
        for line in vf:
            if line.startswith("##"):
                sv_assignment_file.write(line)
                continue
            elif line.startswith("#"):
                sv_assignment_file.write(line)
                sv_assignment_info_file.write(
                    line.rstrip() + "\tNumHp1ReadsFromColumn8\t"
                    "NumHp2ReadsFromColumn8"
                    "\tNumMaternalReadsFromColumn8"
                    "\tNumPaternalReadsFromColumn8\n"
                )
                continue
            line = line.rstrip().split("\t")
            chromosome = line[0]
            format_values = dict(zip(line[8].split(":"), line[9].split(":")))
            if "PS" in format_values:
                del format_values["PS"]

            gt = format_values["GT"]
            replace_patterns = []
            if not include_all_variants and line[6] not in ["PASS", "."]:
                replace_patterns = [("|", "/"), ("1/0", "0/1"), ("2/1", "1/2")]
                additional_values = ["NA"] * 4

            elif chromosome in chrom_hp_origin:
                if gt in ["0/1", "1/0", "0|1", "1|0"] and "RNAMES=" in line[7]:
                    hp_counts = count_reads_haplotypes(
                        chromosome, reads_hap, read_column=line[7]
                    )
                    haplotype = determine_haplotype_from_counts(
                        hp_counts, hapratio, min_read_reassignment
                    )
                    if haplotype == "HP1":
                        if chrom_hp_origin[chromosome]["HP1"][0] == "maternal":
                            format_values["PS"] = "Mat"
                            format_values["GT"] = "1|0"
                            additional_values = [
                                str(hp_counts[1]),
                                str(hp_counts[2]),
                                str(hp_counts[2]),
                                str(hp_counts[1]),
                            ]
                        elif chrom_hp_origin[chromosome]["HP1"][0] == "paternal":
                            format_values["PS"] = "Pat"
                            format_values["GT"] = "0|1"
                            additional_values = [
                                str(hp_counts[1]),
                                str(hp_counts[2]),
                                str(hp_counts[1]),
                                str(hp_counts[2]),
                            ]
                        else:
                            raise RuntimeError("Unknown assignment.")
                    elif haplotype == "HP2":
                        if chrom_hp_origin[chromosome]["HP2"][0] == "maternal":
                            format_values["PS"] = "Mat"
                            format_values["GT"] = "1|0"
                            additional_values = [
                                str(hp_counts[1]),
                                str(hp_counts[2]),
                                str(hp_counts[1]),
                                str(hp_counts[2]),
                            ]

                        elif chrom_hp_origin[chromosome]["HP2"][0] == "paternal":
                            format_values["PS"] = "Pat"
                            format_values["GT"] = "0|1"
                            additional_values = [
                                str(hp_counts[1]),
                                str(hp_counts[2]),
                                str(hp_counts[2]),
                                str(hp_counts[1]),
                            ]
                        else:
                            raise RuntimeError("Unknown assignment.")
                    else:
                        replace_patterns = [("1/0", "0/1")]
                        additional_values = [
                            str(hp_counts[1]),
                            str(hp_counts[2]),
                            "NA",
                            "NA",
                        ]
                else:
                    replace_patterns = [("|", "/"), ("1/0", "0/1"), ("2/1", "1/2")]
                    additional_values = ["NA"] * 4
            else:
                if gt in ["0/1", "1/0", "0|1", "1|0"] and "RNAMES=" in line[7]:
                    hp_counts = count_reads_haplotypes(
                        chromosome, reads_hap, read_column=line[7]
                    )
                    haplotype = determine_haplotype_from_counts(
                        hp_counts, hapratio, min_read_reassignment
                    )

                    additional_values = [
                        str(hp_counts[1]),
                        str(hp_counts[2]),
                        "NA",
                        "NA",
                    ]
                    if haplotype == "HP1":
                        format_values["GT"] = "1|0"
                        format_values["PS"] = haplotype
                    elif haplotype == "HP2":
                        format_values["GT"] = "0|1"
                        format_values["PS"] = haplotype
                    else:
                        replace_patterns = [("1/0", "0/1")]
                else:
                    additional_values = ["NA"] * 4
                    replace_patterns = [("|", "/"), ("1/0", "0/1"), ("2/1", "1/2")]
            write_sv_line(
                line,
                format_values,
                additional_values,
                replace_patterns,
                sv_assignment_file,
                sv_assignment_info_file,
            )


def write_sam_phase(
    in_file: str,
    out_file: str,
    reads_hap: Dict[Tuple[str, str], int],
    mapping_quality: int,
    include_supplementary: bool,
) -> None:
    """Write SAM/BAM file with haplotype tags based on read assignments.

    Args:
        in_file: Input SAM/BAM file path
        out_file: Output file path for gzipped SAM
        reads_hap: Dictionary mapping (chrom, read_id) to haplotype assignment (1 or 2)
        mapping_quality: Minimum mapping quality threshold
        include_supplementary: Whether to include supplementary alignments
    """
    out_nonpofo_bam = gzip.open(out_file, "wb")
    with openfile(in_file) as sf:
        for line in sf:
            if line.startswith(("@HD", "@SQ", "@RG")):
                out_nonpofo_bam.write(line.encode())
                continue
            elif line.startswith(("@")):
                continue
            line = line.rstrip().split("\t")
            out_read = "\t".join(line) + "\n"
            if int(line[4]) < mapping_quality or line[1] not in [
                "0",
                "16",
                "2048",
                "2064",
            ]:
                continue
            elif not include_supplementary and line[1] in ["2048", "2064"]:
                continue
            elif (line[2], line[0]) in reads_hap:
                out_read = (
                    "\t".join(line)
                    + "\t"
                    + "HP:i:"
                    + str(reads_hap[(line[2], line[0])])
                    + "\n"
                )
                out_nonpofo_bam.write(out_read.encode())
    out_nonpofo_bam.close()


def alignment_writer(
    bam: pysam.AlignmentFile,
    chrom: str,
    reads_hap: Dict[Tuple[str, str], int],
    chrom_hp_origin: Dict[str, Dict[str, List[Union[str, int, float]]]],
    outfile: pysam.AlignmentFile,
    mapping_quality: int,
    include_supplementary: bool,
) -> None:
    """Write alignments with correct haplotype tags based on parent of origin.

    Processes alignments for a chromosome and writes them to output with
    haplotype tags set based on parent of origin assignments. Handles maternal/paternal
    assignments differently.

    Args:
        bam: Input BAM file handle
        chrom: Chromosome to process
        reads_hap: Dictionary mapping (chrom, read_id) to haplotype assignment (1 or 2)
        chrom_hp_origin: Dictionary mapping chromosomes to haplotype parent-of-origin assignments.
                        Structure: {chrom: {"HP1": [origin, stats...], "HP2": [origin, stats...]}}
        outfile: Output BAM file handle
        mapping_quality: Minimum mapping quality threshold
        include_supplementary: Whether to include supplementary alignments
    """
    bamiter = bam.fetch(chrom)
    if chrom in chrom_hp_origin and chrom_hp_origin[chrom]["HP1"][0] == "maternal":
        for read in bamiter:
            read.set_tag("HP", None)
            read.set_tag("PS", None)
            if (
                read.mapping_quality < mapping_quality
                or read.is_secondary
                or read.is_qcfail
                or read.is_duplicate
                or read.is_unmapped
                or (read.is_supplementary and not include_supplementary)
            ):
                outfile.write(read)
                continue
            read_id = read.query_name
            ref_name = read.reference_name
            if (ref_name, read_id) in reads_hap:
                read.set_tag("HP", int(reads_hap[(ref_name, read_id)]))
            outfile.write(read)
    elif chrom in chrom_hp_origin and chrom_hp_origin[chrom]["HP2"][0] == "maternal":
        for read in bamiter:
            read.set_tag("HP", None)
            read.set_tag("PS", None)
            if (
                read.mapping_quality < mapping_quality
                or read.is_secondary
                or read.is_qcfail
                or read.is_duplicate
                or read.is_unmapped
                or (read.is_supplementary and not include_supplementary)
            ):
                outfile.write(read)
                continue
            read_id = read.query_name
            ref_name = read.reference_name
            if (ref_name, read_id) in reads_hap:
                if reads_hap[(ref_name, read_id)] == 1:
                    read.set_tag("HP", 2)
                else:
                    read.set_tag("HP", 1)
            outfile.write(read)
    else:
        for read in bamiter:
            read.set_tag("HP", None)
            read.set_tag("PS", None)
            outfile.write(read)
