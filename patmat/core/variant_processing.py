import multiprocessing as mp
import os
import warnings
from collections import defaultdict
from itertools import repeat
import datetime

import pysam
import tabix
from tqdm import tqdm

from patmat.core.phase_processing import get_block
from patmat.io.bam import openalignment
from patmat.io.file_utils import openfile


def process_strand_seq_vcf(
    vcf, chrom, strand_vcf, phased, include_all_variants, ignore_blocks_single
):
    """Process strand-seq VCF file and return phasing information.

    Args:
        args: Command line arguments
        vcf: Path to main VCF file
        chrom: Chromosome being processed

    Returns:
        tuple: (final_dict, strand_phased_vars, phase_block_stat, blocks_dict)
        phase_block_stat and blocks_dict will be None if not phased

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


def get_variant_info(
    feed_list,
    alignment_file,
    chrom,
    mapping_quality,
    include_supplementary,
):
    """
    This function maps each read to heterozygous variants and returns a list of
    haplotype 1, haplotype 2, and unphased variants mapped to each read.
    """
    read_var_list = defaultdict(list)
    samfile = pysam.AlignmentFile(alignment_file, "rb")

    varinfo_dict = {varinfo[1]: varinfo for varinfo in feed_list}
    # lists are sorted, so get min and max position to fetch pileup:
    min_pos = feed_list[0][1]
    max_pos = feed_list[-1][1]
    sam_pileup = samfile.pileup(chrom, min_pos, max_pos + 1, truncate=True)

    for pileupcolumn in sam_pileup:
        # Check if this column is present in our variant list:
        if pileupcolumn.pos not in varinfo_dict:
            continue
        pileupcolumn.set_min_base_quality(0)
        gt, position, ref, alt = varinfo_dict[pileupcolumn.pos]
        for pileupread in pileupcolumn.pileups:
            read_mq = pileupread.alignment.mapping_quality
            if read_mq < mapping_quality:
                continue
            if pileupread.alignment.is_supplementary and not include_supplementary:
                continue
            read_id = pileupread.alignment.query_name
            flag = pileupread.alignment.flag
            ext_key = (read_id, str(flag))  # ,str(read_mq),str(read_start+1))
            read_base = None
            read_base1 = None
            read_base2 = None
            hapt = "NA"
            if pileupread.is_del or pileupread.is_refskip:
                read_var_list[(chrom, position + 1, "noninfo")].append(
                    ":".join((*ext_key, hapt))
                )
                continue
            if gt != "1/2":
                if len(ref) == 1 and len(alt) == 1:  # dealing with mismatches
                    read_base = pileupread.alignment.query_sequence[
                        pileupread.query_position
                    ].upper()
                elif len(ref) > 1 and len(alt) == 1:  # dealing with deletions
                    if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - 1:
                        read_base = pileupread.alignment.query_sequence[
                            pileupread.query_position
                        ].upper()
                    elif pileupread.indel == 0:
                        read_base = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                        ].upper()
                    else:
                        read_var_list[(chrom, position + 1, "noninfo")].append(
                            ":".join((*ext_key, hapt))
                        )
                        continue
                elif len(ref) == 1 and len(alt) > 1:  # dealing with insertions
                    if pileupread.indel > 0 and pileupread.indel == len(alt) - 1:
                        read_base = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(alt)
                        ].upper()
                    elif pileupread.indel == 0:
                        read_base = pileupread.alignment.query_sequence[
                            pileupread.query_position
                        ].upper()
                    else:
                        read_var_list[(chrom, position + 1, "noninfo")].append(
                            ":".join((*ext_key, hapt))
                        )
                        continue

                if read_base == ref:
                    if gt == "1|0":
                        hapt = "2"
                    elif gt == "0|1":
                        hapt = "1"
                elif read_base == alt:
                    if gt == "1|0":
                        hapt = "1"
                    elif gt == "0|1":
                        hapt = "2"
                else:
                    read_base = "noninfo"
                read_var_list[(chrom, position + 1, read_base)].append(
                    ":".join((*ext_key, hapt))
                )

            else:
                alt1, alt2 = alt
                if len(ref) > 1 and len(alt1) > 1 and len(alt2) > 1:
                    read_var_list[(chrom, position + 1, "noninfo")].append(
                        ":".join((*ext_key, hapt))
                    )
                    continue
                if len(ref) == 1 and len(alt1) == 1:  # dealing with mismatches
                    read_base1 = pileupread.alignment.query_sequence[
                        pileupread.query_position
                    ].upper()
                elif len(ref) > 1 and len(ref) > len(alt1):  # dealing with deletions
                    if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - len(
                        alt1
                    ):
                        read_base1 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                            - abs(pileupread.indel)
                        ].upper()
                    elif pileupread.indel == 0:
                        read_base1 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                        ].upper()
                elif len(alt1) > 1 and len(alt1) > len(ref):  # dealing with insertions
                    if pileupread.indel > 0 and pileupread.indel == len(alt1) - len(
                        ref
                    ):
                        read_base1 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(alt1)
                        ].upper()
                    elif pileupread.indel == 0:
                        read_base1 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                        ].upper()
                if len(ref) == 1 and len(alt2) == 1:  # dealing with mismatches
                    read_base2 = pileupread.alignment.query_sequence[
                        pileupread.query_position
                    ].upper()
                elif len(ref) > 1 and len(ref) > len(alt2):  # dealing with deletions
                    if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - len(
                        alt2
                    ):
                        read_base2 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                            - abs(pileupread.indel)
                        ].upper()
                    elif pileupread.indel == 0:
                        read_base2 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                        ].upper()
                elif len(alt2) > 1 and len(alt2) > len(ref):  # dealing with insertions
                    if pileupread.indel > 0 and pileupread.indel == len(alt2) - len(
                        ref
                    ):
                        read_base2 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(alt2)
                        ].upper()
                    elif pileupread.indel == 0:
                        read_base2 = pileupread.alignment.query_sequence[
                            pileupread.query_position : pileupread.query_position
                            + len(ref)
                        ].upper()
                if read_base1 == alt1 or read_base1 == alt2:
                    read_var_list[(chrom, position + 1, read_base1)].append(
                        ":".join((*ext_key, hapt))
                    )
                elif read_base2 == alt1 or read_base2 == alt2:
                    read_var_list[(chrom, position + 1, read_base2)].append(
                        ":".join((*ext_key, hapt))
                    )
                else:
                    read_var_list[(chrom, position + 1, "noninfo")].append(
                        ":".join((*ext_key, hapt))
                    )
    return read_var_list


def write_per_read_variant_file(
    vcf_dict,
    bam_file,
    chunk,
    processes,
    out_file,
    mapping_quality,
    include_supplementary,
):
    """
    This function extracts per-read information for variants.
    """
    outfile = open(out_file, "w")
    for chrom, feed_list in vcf_dict.items():
        # TODO: is this really necessary? only getting the count
        # does it break if bam_file has no alignments? do we see that?
        # bamiter, bam, count = openalignment(bam_file, chrom)
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


def vcf2dict(vcf_strand, chrom):
    """
    Process and converts the strand-seq vcf file to a dictionary
    for downstream use.
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


def strand_vcf2dict_phased(vcf_strand, vcf, include_all_variants, chrom):
    """
    Intersects input vcf and strand-seq vcf and stores phased
    and unphased heterozygous variants into a dictionary for read phasing.
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
    blocks, vcf_strand, vcf, chrom, include_all_variants, ignore_blocks_single
):
    """
    In case --phased option is given, this function uses phased variants
    from strand-seq to correct PofO phasing switches across phased blocks
    . This function then intersects input vcf and strand-seq vcf and
    stores phased and unphased heterozygous variants into a
    dictionary for read phasing.
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
