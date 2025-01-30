import gzip
import os
import warnings
from collections import defaultdict
from typing import DefaultDict, Dict, List, Set, Tuple

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


def pofo_sv_write(
    sv_file,
    out,
    chrom_hp_origin,
    reads_hap,
    min_read_reassignment,
    include_all_variants,
    hapratio,
):
    sv_assignment_file = open(
        out + "_" + os.path.basename(sv_file) + "_PofO_Assignment_SVs.vcf", "w"
    )
    sv_assignment_file_info = open(
        out + "_" + os.path.basename(sv_file) + "_Variant_Assignment_SVs_info.tsv", "w"
    )
    with openfile(sv_file) as vf:
        for line in vf:
            if line.startswith("##"):
                sv_assignment_file.write(line)
                continue
            elif line.startswith("#"):
                sv_assignment_file.write(line)
                sv_assignment_file_info.write(
                    line.rstrip() + "\tNumHp1ReadsFromColumn8\t"
                    "NumHp2ReadsFromColumn8"
                    "\tNumMaternalReadsFromColumn8"
                    "\tNumPaternalReadsFromColumn8\n"
                )
                continue
            line = line.rstrip().split("\t")
            gt = line[9].split(":")[0]
            if "PS" in line[8].split(":"):
                ps_index = line[8].split(":").index("PS")
                new_ps = line[8].split(":")
                new_hp = line[9].split(":")
                new_ps.pop(ps_index)
                new_hp.pop(ps_index)
            else:
                new_ps = line[8].split(":")
                new_hp = line[9].split(":")
            if not include_all_variants and line[6] not in ["PASS", "."]:
                line_out = (
                    line[0:8]
                    + [":".join(new_ps)]
                    + [
                        ":".join(new_hp)
                        .replace("|", "/")
                        .replace("1/0", "0/1")
                        .replace("2/1", "1/2")
                    ]
                    + ["NA"] * 4
                )
                sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                sv_assignment_file_info.write("\t".join(line_out) + "\n")
                continue
            if line[0] in chrom_hp_origin:
                if gt in ["0/1", "1/0", "0|1", "1|0"] and "RNAMES=" in line[7]:
                    hp1_count = 0
                    hp2_count = 0
                    for read_ID in line[7].split("RNAMES=")[1].split(";")[0].split(","):
                        if (line[0], read_ID) in reads_hap:
                            if reads_hap[(line[0], read_ID)] == 1:
                                hp1_count += 1
                            elif reads_hap[(line[0], read_ID)] == 2:
                                hp2_count += 1
                    if (
                        hp1_count > hp2_count
                        and hp1_count / (hp1_count + hp2_count) >= hapratio
                        and hp1_count >= min_read_reassignment
                    ):
                        if chrom_hp_origin[line[0]]["HP1"][0] == "maternal":
                            line_out = (
                                line[0:8]
                                + [":".join(new_ps) + ":PS"]
                                + ["1|0:" + ":".join(new_hp[1:]) + ":Mat"]
                                + [
                                    str(hp1_count),
                                    str(hp2_count),
                                    str(hp2_count),
                                    str(hp1_count),
                                ]
                            )
                            sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                            sv_assignment_file_info.write("\t".join(line_out) + "\n")
                        if chrom_hp_origin[line[0]]["HP1"][0] == "paternal":
                            line_out = (
                                line[0:8]
                                + [":".join(new_ps) + ":PS"]
                                + ["0|1:" + ":".join(new_hp[1:]) + ":Pat"]
                                + [
                                    str(hp1_count),
                                    str(hp2_count),
                                    str(hp1_count),
                                    str(hp2_count),
                                ]
                            )
                            sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                            sv_assignment_file_info.write("\t".join(line_out) + "\n")
                    elif (
                        hp2_count > hp1_count
                        and hp2_count / (hp1_count + hp2_count) >= hapratio
                        and hp2_count >= min_read_reassignment
                    ):
                        if chrom_hp_origin[line[0]]["HP2"][0] == "maternal":
                            line_out = (
                                line[0:8]
                                + [":".join(new_ps) + ":PS"]
                                + ["1|0:" + ":".join(new_hp[1:]) + ":Mat"]
                                + [
                                    str(hp1_count),
                                    str(hp2_count),
                                    str(hp1_count),
                                    str(hp2_count),
                                ]
                            )
                            sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                            sv_assignment_file_info.write("\t".join(line_out) + "\n")
                        if chrom_hp_origin[line[0]]["HP2"][0] == "paternal":
                            line_out = (
                                line[0:8]
                                + [":".join(new_ps) + ":PS"]
                                + ["0|1:" + ":".join(new_hp[1:]) + ":Pat"]
                                + [
                                    str(hp1_count),
                                    str(hp2_count),
                                    str(hp2_count),
                                    str(hp1_count),
                                ]
                            )
                            sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                            sv_assignment_file_info.write("\t".join(line_out) + "\n")
                    else:
                        line_out = (
                            line[0:8]
                            + [":".join(new_ps)]
                            + [line[9].replace("|", "/").replace("1/0", "0/1")]
                            + [str(hp1_count), str(hp2_count), "NA", "NA"]
                        )
                        sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                        sv_assignment_file_info.write("\t".join(line_out) + "\n")
                else:
                    line_out = (
                        line[0:8]
                        + [":".join(new_ps)]
                        + [
                            line[9]
                            .replace("|", "/")
                            .replace("1/0", "0/1")
                            .replace("2/1", "1/2")
                        ]
                        + ["NA"] * 4
                    )
                    sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                    sv_assignment_file_info.write("\t".join(line_out) + "\n")
            else:
                if gt in ["0/1", "1/0", "0|1", "1|0"] and "RNAMES=" in line[7]:
                    hp1_count = 0
                    hp2_count = 0
                    for read_ID in line[7].split("RNAMES=")[1].split(";")[0].split(","):
                        if (line[0], read_ID) in reads_hap:
                            if reads_hap[(line[0], read_ID)] == 1:
                                hp1_count += 1
                            elif reads_hap[(line[0], read_ID)] == 2:
                                hp2_count += 1
                    if (
                        hp1_count > hp2_count
                        and hp1_count / (hp1_count + hp2_count) >= hapratio
                        and hp1_count >= min_read_reassignment
                    ):
                        line_out = (
                            line[0:8]
                            + [":".join(new_ps) + ":PS"]
                            + ["1|0:" + ":".join(new_hp[1:]) + ":HP1"]
                            + [str(hp1_count), str(hp2_count), "NA", "NA"]
                        )
                        sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                        sv_assignment_file_info.write("\t".join(line_out) + "\n")
                    elif (
                        hp2_count > hp1_count
                        and hp2_count / (hp1_count + hp2_count) >= hapratio
                        and hp2_count >= min_read_reassignment
                    ):
                        line_out = (
                            line[0:8]
                            + [":".join(new_ps) + ":PS"]
                            + ["0|1:" + ":".join(new_hp[1:]) + ":HP2"]
                            + [str(hp1_count), str(hp2_count), "NA", "NA"]
                        )
                        sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                        sv_assignment_file_info.write("\t".join(line_out) + "\n")
                    else:
                        line_out = (
                            line[0:8]
                            + [":".join(new_ps)]
                            + [line[9].replace("|", "/").replace("1/0", "0/1")]
                            + [str(hp1_count), str(hp2_count), "NA", "NA"]
                        )
                        sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                        sv_assignment_file_info.write("\t".join(line_out) + "\n")
                else:
                    line_out = (
                        line[0:8]
                        + [":".join(new_ps)]
                        + [
                            line[9]
                            .replace("|", "/")
                            .replace("1/0", "0/1")
                            .replace("2/1", "1/2")
                        ]
                        + ["NA"] * 4
                    )
                    sv_assignment_file.write("\t".join(line_out[0:-4]) + "\n")
                    sv_assignment_file_info.write("\t".join(line_out) + "\n")


def write_sam_phase(
    in_file, out_file, reads_hap, mapping_quality, include_supplementary
):
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
    bam,
    chrom,
    reads_hap,
    chrom_hp_origin,
    outfile,
    mapping_quality,
    include_supplementary,
):
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
