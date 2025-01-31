import multiprocessing as mp
import os
import warnings
from collections import defaultdict
from itertools import repeat
from typing import Any, DefaultDict, Dict, List, Optional, Sequence, Set, Tuple, Union

import pysam
import tabix
from tqdm import tqdm

from patmat.core.phase_processing import get_block
from patmat.io.file_utils import openfile


def process_strand_seq_vcf(
    vcf: str,
    chrom: str,
    strand_vcf: Optional[str],
    phased: bool,
    include_all_variants: bool,
    ignore_blocks_single: bool,
) -> Tuple[
    DefaultDict[str, Set[Tuple]],
    int,
    Optional[DefaultDict[Tuple[str, str, str], int]],
    Optional[Dict[Tuple[str, str], Tuple[str, str]]],
]:
    """Process strand-seq VCF file and return phasing information.

    Args:
        vcf: Path to main VCF file
        chrom: Chromosome being processed
        strand_vcf: Path to strand-seq VCF file
        phased: Whether input VCF is phased
        include_all_variants: Whether to include variants not marked as PASS
        ignore_blocks_single: Whether to ignore single variant blocks

    Returns:
        Tuple containing:
            - Dictionary mapping chromosomes to sets of variant tuples
            - Number of strand phased variants
            - Dictionary of phase block statistics (None if not phased)
            - Dictionary of block boundaries (None if not phased)

    Raises:
        Exception: If no strand-seq VCF is provided
    """
    if strand_vcf is None:
        raise Exception("No strand-seq VCF is given.")

    vcf_strand = os.path.abspath(strand_vcf)

    if not phased:
        final_dict, strand_phased_vars = strand_vcf2dict_phased(
            vcf_strand, vcf, include_all_variants, chrom
        )
        return final_dict, strand_phased_vars, None, None

    else:
        blocks_phased, blocks_dict = get_block(vcf, chrom)
        final_dict, strand_phased_vars, phase_block_stat = vcf2dict_phased(
            blocks_phased,
            vcf_strand,
            vcf,
            chrom,
            include_all_variants,
            ignore_blocks_single,
        )
        return final_dict, strand_phased_vars, phase_block_stat, blocks_dict


def process_read_base(
    pileupread: pysam.PileupRead, ref: str, alt: str, query_position: int
) -> Optional[str]:
    """Process and return read base for single alternative.

    Args:
        pileupread: Pileup read object
        ref: Reference allele
        alt: Alternative allele
        query_position: Position in query sequence

    Returns:
        Processed base string or None if invalid/unprocessable
    """
    if len(ref) == 1 and len(alt) == 1:  # mismatch
        return pileupread.alignment.query_sequence[query_position].upper()

    elif len(ref) > 1 and len(alt) == 1:  # deletion
        if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - 1:
            return pileupread.alignment.query_sequence[query_position].upper()
        elif pileupread.indel == 0:
            return pileupread.alignment.query_sequence[
                query_position : query_position + len(ref)
            ].upper()

    elif len(ref) == 1 and len(alt) > 1:  # insertion
        if pileupread.indel > 0 and pileupread.indel == len(alt) - 1:
            return pileupread.alignment.query_sequence[
                query_position : query_position + len(alt)
            ].upper()
        elif pileupread.indel == 0:
            return pileupread.alignment.query_sequence[query_position].upper()

    return None


def process_complex_read_base(
    pileupread: pysam.PileupRead, ref: str, alt: str, query_position: int
) -> Optional[str]:
    """Process and return read base for complex variants.

    Args:
        pileupread: Pileup read object
        ref: Reference allele
        alt: Alternative allele
        query_position: Position in query sequence

    Returns:
        Processed base string or None if invalid/unprocessable
    """
    if len(ref) == 1 and len(alt) == 1:  # mismatch
        return pileupread.alignment.query_sequence[query_position].upper()

    elif len(ref) > 1 and len(ref) > len(alt):  # deletion
        if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - len(alt):
            return pileupread.alignment.query_sequence[
                query_position : query_position + len(ref) - abs(pileupread.indel)
            ].upper()
        elif pileupread.indel == 0:
            return pileupread.alignment.query_sequence[
                query_position : query_position + len(ref)
            ].upper()

    elif len(alt) > 1 and len(alt) > len(ref):  # insertion
        if pileupread.indel > 0 and pileupread.indel == len(alt) - len(ref):
            return pileupread.alignment.query_sequence[
                query_position : query_position + len(alt)
            ].upper()
        elif pileupread.indel == 0:
            return pileupread.alignment.query_sequence[
                query_position : query_position + len(ref)
            ].upper()

    return None


def process_simple_variant(
    pileupread: pysam.PileupRead,
    ref: str,
    alt: str,
    gt: str,
    chrom: str,
    position: int,
    ext_key: Tuple[str, str],
) -> Tuple[str, str, Tuple[str, int, str]]:
    """Process a simple variant (not 1/2) and return base, haplotype, and key.

    Args:
        pileupread: Pileup read object
        ref: Reference allele
        alt: Alternative allele
        gt: Genotype string
        chrom: Chromosome name
        position: Variant position
        ext_key: Tuple of (read_name, flag)

    Returns:
        Tuple containing:
            - Read base string
            - Haplotype assignment ("1", "2", or "NA")
            - Key tuple (chrom, position+1, read_base)
    """
    read_base = process_read_base(pileupread, ref, alt, pileupread.query_position)
    if read_base is None:
        return "noninfo", "NA", (chrom, position + 1, "noninfo")

    hapt = "NA"
    if read_base == ref:
        hapt = "2" if gt == "1|0" else "1" if gt == "0|1" else "NA"
    elif read_base == alt:
        hapt = "1" if gt == "1|0" else "2" if gt == "0|1" else "NA"
    else:
        read_base = "noninfo"

    return read_base, hapt, (chrom, position + 1, read_base)


def process_complex_variant(
    pileupread: pysam.PileupRead,
    ref: str,
    alt: Tuple[str, str],
    chrom: str,
    position: int,
    ext_key: Tuple[str, str],
) -> Tuple[str, str, Tuple[str, int, str]]:
    """Process a complex variant (1/2) and return base, haplotype, and key.

    Args:
        pileupread: Pileup read object
        ref: Reference allele
        alt: Tuple of two alternative alleles
        chrom: Chromosome name
        position: Variant position
        ext_key: Tuple of (read_name, flag)

    Returns:
        Tuple containing:
            - Read base string
            - Haplotype assignment ("NA" for complex variants)
            - Key tuple (chrom, position+1, read_base)
    """
    alt1, alt2 = alt
    if len(ref) > 1 and len(alt1) > 1 and len(alt2) > 1:
        return "noninfo", "NA", (chrom, position + 1, "noninfo")

    read_base1 = process_complex_read_base(
        pileupread, ref, alt1, pileupread.query_position
    )
    read_base2 = process_complex_read_base(
        pileupread, ref, alt2, pileupread.query_position
    )

    if read_base1 and (read_base1 == alt1 or read_base1 == alt2):
        return read_base1, "NA", (chrom, position + 1, read_base1)
    elif read_base2 and (read_base2 == alt1 or read_base2 == alt2):
        return read_base2, "NA", (chrom, position + 1, read_base2)
    else:
        return "noninfo", "NA", (chrom, position + 1, "noninfo")


def get_variant_info(
    feed_list: Sequence[Tuple[str, int, str, Union[str, Tuple[str, str]]]],
    alignment_file: str,
    chrom: str,
    mapping_quality: int,
    include_supplementary: bool,
) -> DefaultDict[Tuple[str, int, str], Sequence[str]]:
    """Map each read to heterozygous variants and return haplotype assignments.

    Args:
        feed_list: List of variant information tuples containing:
            (genotype, position, reference, alternative)
        alignment_file: Path to the alignment file
        chrom: Chromosome name
        mapping_quality: Minimum mapping quality threshold
        include_supplementary: Whether to include supplementary alignments

    Returns:
        Dictionary mapping (chrom, pos, base) to sequence of read information strings
    """
    read_var_list = defaultdict(list)
    samfile = pysam.AlignmentFile(alignment_file, "rb")

    varinfo_dict = {varinfo[1]: varinfo for varinfo in feed_list}
    min_pos, max_pos = feed_list[0][1], feed_list[-1][1]

    for pileupcolumn in samfile.pileup(chrom, min_pos, max_pos + 1, truncate=True):
        if pileupcolumn.pos not in varinfo_dict:
            continue

        pileupcolumn.set_min_base_quality(0)
        gt, position, ref, alt = varinfo_dict[pileupcolumn.pos]

        for pileupread in pileupcolumn.pileups:
            # Filter reads based on quality and supplementary status
            if pileupread.alignment.mapping_quality < mapping_quality:
                continue
            if pileupread.alignment.is_supplementary and not include_supplementary:
                continue

            # Create extended key for read identification
            ext_key = (pileupread.alignment.query_name, str(pileupread.alignment.flag))

            # Handle deletions and reference skips
            if pileupread.is_del or pileupread.is_refskip:
                read_var_list[(chrom, position + 1, "noninfo")].append(
                    ":".join((*ext_key, "NA"))
                )
                continue

            # Process variant based on genotype
            if gt != "1/2":
                read_base, hapt, key = process_simple_variant(
                    pileupread, ref, alt, gt, chrom, position, ext_key
                )
            else:
                read_base, hapt, key = process_complex_variant(
                    pileupread, ref, alt, chrom, position, ext_key
                )

            read_var_list[key].append(":".join((*ext_key, hapt)))

    return read_var_list


def write_per_read_variant_file(
    vcf_dict: Dict[str, Sequence[Tuple[Any, ...]]],
    bam_file: str,
    chunk: int,
    processes: int,
    out_file: str,
    mapping_quality: int,
    include_supplementary: bool,
) -> None:
    """Extract and write per-read information for variants.

    Args:
        vcf_dict: Dictionary mapping chromosomes to sequences of variant tuples
        bam_file: Path to BAM/CRAM file
        chunk: Size of chunks for parallel processing
        processes: Number of parallel processes
        out_file: Output file path
        mapping_quality: Minimum mapping quality threshold
        include_supplementary: Whether to include supplementary alignments
    """
    outfile = open(out_file, "w")
    for chrom, feed_list in vcf_dict.items():
        # TODO: is this really necessary? only getting the count
        # does it break if bam_file has no alignments? do we see that?
        # bamiter, bam, count = open_alignment(bam_file, chrom)
        if True:  #  count > 0:
            # sort by position:
            feed_list = sorted(feed_list, key=lambda x: x[1])
            vcf_info_list = [
                list(feed_list)[x : x + chunk] for x in range(0, len(feed_list), chunk)
            ]

            p = mp.Pool(processes)
            results = p.starmap(
                get_variant_info,
                list(
                    zip(
                        vcf_info_list,
                        repeat(bam_file),
                        repeat(chrom),
                        repeat(mapping_quality),
                        repeat(include_supplementary),
                    )
                ),
            )
            p.close()
            p.join()

            for result in results:
                if result is not None:
                    for key, val in result.items():
                        outfile.write(
                            "\t".join(map(str, key)) + "\t" + ",".join(val) + "\n"
                        )

        else:
            warnings.warn(
                "{} does not have any mapped reads in alignment "
                "file Or alignment is truncated or corrupt indexed. "
                "Skipping it.".format(chrom)
            )

    outfile.close()


def vcf2dict(vcf_strand: str, chrom: str) -> DefaultDict[str, Dict[str, str]]:
    """Process and convert the strand-seq VCF file to a dictionary.

    Args:
        vcf_strand: Path to strand-seq VCF file
        chrom: Chromosome to process

    Returns:
        Nested defaultdict mapping chromosome to position to genotype string

    Notes:
        Handles special cases for partially phased genotypes (.|0, 0|., etc)
        converting them to fully phased representation
    """
    vcf_dict = defaultdict(dict)
    vcf_file = openfile(vcf_strand)
    for line in vcf_file:
        line = line.rstrip().split("\t")
        if line[0] != chrom:
            continue
        if line[9].startswith((".|0", "0|.", "1|.", ".|1")):
            line[9] = (
                line[9]
                .replace(".|0", "1|0")
                .replace(".|1", "0|1")
                .replace("1|.", "1|0")
                .replace("0|.", "0|1")
            )
            warnings.warn(
                "{}:{} variant in strand-seq vcf has .|0 or 0|. "
                "or .|1 or 1|. phased genotype. Note that it "
                "will be interpreted as 1|0 or 0|1 or 0|1 or "
                "1|0".format(line[0], line[1])
            )
        if line[9].startswith(("1|0", "0|1")):
            vcf_dict[line[0]][line[1]] = line[9].split(":")[0]
    vcf_file.close()
    return vcf_dict


def strand_vcf2dict_phased(
    vcf_strand: str, vcf: str, include_all_variants: bool, chrom: str
) -> Tuple[
    DefaultDict[str, Set[Tuple[str, int, str, Union[str, Tuple[str, str]]]]], int
]:
    """Intersect input VCF and strand-seq VCF to store variants for read phasing.

    Args:
        vcf_strand: Path to strand-seq VCF file
        vcf: Path to input VCF file
        include_all_variants: Whether to include variants not marked as PASS
        chrom: Chromosome to process

    Returns:
        Tuple containing:
            - Dictionary mapping chromosomes to sets of variant tuples
            - Number of strand phased variants
    """
    final_dict = defaultdict(set)
    vcf_dict = vcf2dict(vcf_strand, chrom)
    with openfile(vcf) as vf:
        for line in vf:
            line = line.rstrip().split("\t")
            if line[0] != chrom:
                continue
            if not include_all_variants and line[6] not in ["PASS", "."]:
                continue
            if line[9].startswith(("0/1", "1/0", "0|1", "1|0")) and (
                line[0] in vcf_dict and line[1] in vcf_dict[line[0]]
            ):
                final_dict[line[0]].add(
                    (
                        vcf_dict[line[0]][line[1]],
                        int(line[1]) - 1,
                        line[3].upper(),
                        line[4].upper(),
                    )
                )
            elif line[9].startswith(("0/1", "1/0", "0|1", "1|0")):
                final_dict[line[0]].add(
                    ("0/1", int(line[1]) - 1, line[3].upper(), line[4].upper())
                )
            elif line[9].startswith(("1/2", "1|2", "2/1", "2|1")):
                final_dict[line[0]].add(
                    (
                        "1/2",
                        int(line[1]) - 1,
                        line[3].upper(),
                        (line[4].split(",")[0].upper(), line[4].split(",")[1].upper()),
                    )
                )
    strand_phased_vars = len(vcf_dict.keys())
    vcf_dict.clear()
    return final_dict, strand_phased_vars


def vcf2dict_phased(
    blocks: Sequence[Tuple[str, int, int]],
    vcf_strand: str,
    vcf: str,
    chrom: str,
    include_all_variants: bool,
    ignore_blocks_single: bool,
) -> Tuple[
    DefaultDict[str, Set[Tuple[str, int, str, Union[str, Tuple[str, str]]]]],
    int,
    DefaultDict[Tuple[str, str, str], int],
]:
    """Process phased variants using strand-seq data for switch correction.

    Process phased blocks to correct switches using strand-seq data, then
    intersect with input VCF for read phasing.

    Args:
        blocks: Sequence of (chrom, start, end) tuples defining phase blocks
        vcf_strand: Path to strand-seq VCF file
        vcf: Path to input VCF file
        chrom: Chromosome to process
        include_all_variants: Whether to include variants not marked as PASS
        ignore_blocks_single: Whether to ignore single variant blocks

    Returns:
        Tuple containing:
            - Dictionary mapping chromosomes to sets of variant tuples
            - Number of strand phased variants
            - Dictionary of phase block statistics tracking agreement/disagreement counts

    Notes:
        Uses sliding window of two variants to detect and correct switches.
        Handles both single variant blocks and blocks with multiple variants.
    """
    final_dict = defaultdict(set)
    phase_block_stat = defaultdict(int)
    vcf_dict = vcf2dict(vcf_strand, chrom)
    tb_vcf = tabix.open(os.path.abspath(vcf))
    for block in blocks:
        b_chrom, b_start, b_end = block
        if b_chrom != chrom:
            continue
        try:
            records_whats = tb_vcf.query(b_chrom, b_start - 1, b_end + 1)
            records_whats = list(records_whats)
        except:
            warnings.warn(
                "{}:{}-{} block cannot be extracted from vcf file. "
                "Make sure file is indexed. Skipping it.".format(
                    b_chrom, b_start, b_end
                )
            )
            continue
        if b_chrom in vcf_dict:
            agreement_count_temp = 0
            disagreement_count_temp = 0
            agreement_count = 0
            disagreement_count = 0
            vars_list = list()
            two_or_more = False
            for vcf_line in records_whats:
                if not include_all_variants and vcf_line[6] not in ["PASS", "."]:
                    continue
                if agreement_count_temp + disagreement_count_temp == 2:
                    two_or_more = True
                    if agreement_count_temp == 2:
                        for var in vars_list:
                            vcf_dict[var[0]][var[1]] = var[9].split(":")[0]
                    elif disagreement_count_temp == 2:
                        for var in vars_list:
                            vcf_dict[var[0]][var[1]] = var[9].split(":")[0][::-1]
                    agreement_count_temp = 0
                    disagreement_count_temp = 0
                    vars_list = list()
                if vcf_line[9].startswith(("1|0", "0|1")):
                    vars_list.append(vcf_line)
                if vcf_line[1] in vcf_dict[vcf_line[0]]:
                    if (
                        vcf_line[9].startswith("1|0")
                        and vcf_dict[vcf_line[0]][vcf_line[1]] == "1|0"
                    ):
                        agreement_count += 1
                        agreement_count_temp += 1
                    elif (
                        vcf_line[9].startswith("1|0")
                        and vcf_dict[vcf_line[0]][vcf_line[1]] == "0|1"
                    ):
                        disagreement_count += 1
                        disagreement_count_temp += 1
                    elif (
                        vcf_line[9].startswith("0|1")
                        and vcf_dict[vcf_line[0]][vcf_line[1]] == "0|1"
                    ):
                        agreement_count += 1
                        agreement_count_temp += 1
                    elif (
                        vcf_line[9].startswith("0|1")
                        and vcf_dict[vcf_line[0]][vcf_line[1]] == "1|0"
                    ):
                        disagreement_count += 1
                        disagreement_count_temp += 1
            if agreement_count_temp + disagreement_count_temp == 1 and (
                two_or_more or not ignore_blocks_single
            ):
                if agreement_count_temp == 1:
                    for var in vars_list:
                        vcf_dict[var[0]][var[1]] = var[9].split(":")[0]
                elif disagreement_count_temp == 1:
                    for var in vars_list:
                        vcf_dict[var[0]][var[1]] = var[9].split(":")[0][::-1]

            if agreement_count > disagreement_count:
                phase_block_stat[
                    (b_chrom, str(b_start), "agreement")
                ] += agreement_count
                phase_block_stat[
                    (b_chrom, str(b_start), "disagreement")
                ] += disagreement_count
            else:
                phase_block_stat[
                    (b_chrom, str(b_start), "agreement")
                ] += disagreement_count
                phase_block_stat[
                    (b_chrom, str(b_start), "disagreement")
                ] += agreement_count

    with openfile(vcf) as vf:
        for line in vf:
            line = line.rstrip().split("\t")
            if line[0] != chrom:
                continue
            if not include_all_variants and line[6] not in ["PASS", "."]:
                continue
            if line[9].startswith(("0/1", "1/0", "0|1", "1|0")) and (
                line[0] in vcf_dict and line[1] in vcf_dict[line[0]]
            ):
                final_dict[line[0]].add(
                    (
                        vcf_dict[line[0]][line[1]],
                        int(line[1]) - 1,
                        line[3].upper(),
                        line[4].upper(),
                    )
                )
            elif line[9].startswith(("0/1", "1/0", "0|1", "1|0")):
                final_dict[line[0]].add(
                    ("0/1", int(line[1]) - 1, line[3].upper(), line[4].upper())
                )
            elif line[9].startswith(("1/2", "1|2", "2/1", "2|1")):
                final_dict[line[0]].add(
                    (
                        "1/2",
                        int(line[1]) - 1,
                        line[3].upper(),
                        (line[4].split(",")[0].upper(), line[4].split(",")[1].upper()),
                    )
                )

    strand_phased_vars = len(vcf_dict.keys())
    vcf_dict.clear()
    return final_dict, strand_phased_vars, phase_block_stat
