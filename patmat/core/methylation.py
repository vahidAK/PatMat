import os
import warnings
from collections import defaultdict
from typing import Any, DefaultDict, Dict, List, TextIO, Tuple, Union

import tabix

from patmat.core.subprocess_calls import run_command
from patmat.io.file_utils import openfile


def init_output_files(out: str, known_dmr: str) -> TextIO:
    """Initialize output files and write headers.

    Args:
        out: Output file prefix
        known_dmr: Path to DMR file

    Returns:
        File handle for methylation output file
    """
    out_meth = open(out + "_CpG-Methylation-Status-at-DMRs.tsv", "w")
    dmr_file = openfile(known_dmr)
    header = next(dmr_file).rstrip()
    dmr_file.close()

    out_meth.write(
        header + "\tAll_CpGs_At_iDMR_CouldBeExaminedInBothHaplotypes\t"
        "DifferentiallyMethylatedCpGs_HypermethylatedOnHP1\t"
        "DifferentiallyMethylatedCpGs_HypermethylatedOnHP2\t"
        "MeanDiffMethylationOf_DifferentiallyMethylatedCpGs_HypermethylatedOnHP1\t"
        "MeanDiffMethylationOf_DifferentiallyMethylatedCpGs_HypermethylatedOnHP2\t"
        "MethylationFrequency_HP1\tMethylationFrequency_HP2\t"
        "Included_Or_Ignored_For_PofO_Assignment\n"
    )
    return out_meth


def check_dmr_validity(
    num_cg: int, diff_counts: Tuple[int, int], min_cg: int, cpg_difference: float
) -> str:
    """Check if DMR meets validity criteria.

    Args:
        num_cg: Number of CpGs in DMR
        diff_counts: Tuple of (diff_cg_hp1, diff_cg_hp2)
        min_cg: Minimum required CpGs
        cpg_difference: Minimum required methylation difference

    Returns:
        Status string indicating validity and reason if invalid
    """
    diff_cg_hp1, diff_cg_hp2 = diff_counts

    if num_cg < 1:
        return "Ignored:No_CpG"

    if num_cg < min_cg and abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
        return "Ignored:DidNotMeet_min_cg_And_cpg_difference"

    if num_cg < min_cg and diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
        return "Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes"

    if num_cg < min_cg:
        return "Ignored:DidNotMeet_min_cg"

    if abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
        return "Ignored:DidNotMeet_cpg_difference"

    if diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
        return "Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes"

    return "Included"


def calculate_diff_scores(
    meth_stats: Tuple[int, int, float, float]
) -> Tuple[float, float, int, int]:
    """Calculate differential methylation scores.

    Args:
        meth_stats: Tuple containing (diff_cg_hp1, diff_cg_hp2, diff_cg_hp1_meth, diff_cg_hp2_meth)

    Returns:
        Tuple of (hp1_score, hp2_score, diff_cg_hp1, diff_cg_hp2)
    """
    diff_cg_hp1, diff_cg_hp2, diff_cg_hp1_meth, diff_cg_hp2_meth = meth_stats
    return (
        diff_cg_hp1 * diff_cg_hp1_meth,
        diff_cg_hp2 * diff_cg_hp2_meth,
        diff_cg_hp1,
        diff_cg_hp2,
    )


def update_maternal_counts(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    num_cg: int,
    diff_scores: Tuple[float, float, int, int],
) -> None:
    """Update counts for maternal DMRs.

    Args:
        chrom_hp_origin_count: Dictionary to update
        dmr_chrom: Chromosome
        num_cg: Number of CpGs
        diff_scores: Tuple of (hp1_score, hp2_score, diff_cg_hp1, diff_cg_hp2)
    """
    hp1_score, hp2_score, diff_cg_hp1, diff_cg_hp2 = diff_scores

    # Update scores
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["score_HP1"] += hp1_score
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["score_HP2"] += hp2_score
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["score_HP2"] += hp1_score
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["score_HP1"] += hp2_score

    # Update differences
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP1"] += diff_cg_hp1
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP2"] += diff_cg_hp2
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP2"] += diff_cg_hp1
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP1"] += diff_cg_hp2

    # Update DMR counts based on scores
    if hp1_score > hp2_score:
        update_maternal_hp1_dominant(
            chrom_hp_origin_count, dmr_chrom, num_cg, diff_cg_hp1, diff_cg_hp2
        )
    elif hp1_score < hp2_score:
        update_maternal_hp2_dominant(
            chrom_hp_origin_count, dmr_chrom, num_cg, diff_cg_hp1, diff_cg_hp2
        )


def update_paternal_counts(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    num_cg: int,
    diff_scores: Tuple[float, float, int, int],
) -> None:
    """Update counts for paternal DMRs.

    Args:
        chrom_hp_origin_count: Dictionary to update
        dmr_chrom: Chromosome
        num_cg: Number of CpGs
        diff_scores: Tuple of (hp1_score, hp2_score, diff_cg_hp1, diff_cg_hp2)
    """
    hp1_score, hp2_score, diff_cg_hp1, diff_cg_hp2 = diff_scores

    # Update scores
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["score_HP1"] += hp2_score
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["score_HP2"] += hp1_score
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["score_HP2"] += hp2_score
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["score_HP1"] += hp1_score

    # Update differences
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP1"] += diff_cg_hp2
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP2"] += diff_cg_hp1
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP2"] += diff_cg_hp2
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP1"] += diff_cg_hp1

    # Update DMR counts based on scores
    if hp1_score > hp2_score:
        update_paternal_hp1_dominant(
            chrom_hp_origin_count, dmr_chrom, num_cg, diff_cg_hp1, diff_cg_hp2
        )
    elif hp1_score < hp2_score:
        update_paternal_hp2_dominant(
            chrom_hp_origin_count, dmr_chrom, num_cg, diff_cg_hp1, diff_cg_hp2
        )


def update_maternal_hp1_dominant(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    num_cg: int,
    diff_cg_hp1: int,
    diff_cg_hp2: int,
) -> None:
    """Update counts when HP1 is dominant in maternal DMR."""
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["dmr_HP1"] += 1
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["dmr_HP2"] += 1
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["allcg_HP1"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["allcg_HP2"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["alldiff_HP1"] += (
        diff_cg_hp1 + diff_cg_hp2
    )
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["alldiff_HP2"] += (
        diff_cg_hp1 + diff_cg_hp2
    )


def update_maternal_hp2_dominant(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    num_cg: int,
    diff_cg_hp1: int,
    diff_cg_hp2: int,
) -> None:
    """Update counts when HP2 is dominant in maternal DMR."""
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["dmr_HP2"] += 1
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["dmr_HP1"] += 1
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["allcg_HP2"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["allcg_HP1"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["alldiff_HP2"] += (
        diff_cg_hp1 + diff_cg_hp2
    )
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["alldiff_HP1"] += (
        diff_cg_hp1 + diff_cg_hp2
    )


def update_paternal_hp1_dominant(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    num_cg: int,
    diff_cg_hp1: int,
    diff_cg_hp2: int,
) -> None:
    """Update counts when HP1 is dominant in paternal DMR."""
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["dmr_HP2"] += 1
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["dmr_HP1"] += 1
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["allcg_HP2"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["allcg_HP1"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["alldiff_HP2"] += (
        diff_cg_hp1 + diff_cg_hp2
    )
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["alldiff_HP1"] += (
        diff_cg_hp1 + diff_cg_hp2
    )


def update_paternal_hp2_dominant(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    num_cg: int,
    diff_cg_hp1: int,
    diff_cg_hp2: int,
) -> None:
    """Update counts when HP2 is dominant in paternal DMR."""
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["dmr_HP1"] += 1
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["dmr_HP2"] += 1
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["allcg_HP1"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["allcg_HP2"] += num_cg
    chrom_hp_origin_count[(dmr_chrom, "maternal")]["alldiff_HP1"] += (
        diff_cg_hp1 + diff_cg_hp2
    )
    chrom_hp_origin_count[(dmr_chrom, "paternal")]["alldiff_HP2"] += (
        diff_cg_hp1 + diff_cg_hp2
    )


def process_dmr_line(line: List[str]) -> Tuple[str, int, int, str]:
    """Process a single DMR line and extract key fields.

    Args:
        line: Split line from DMR file

    Returns:
        Tuple of (chromosome, start position, end position, origin)
    """
    dmr_chrom = line[0]
    dmr_start = int(line[1]) - 1
    dmr_end = int(line[2])
    origin = line[3].lower()

    if origin not in ["maternal", "paternal"]:
        idmr = "\t".join(line)
        warnings.warn(
            f"iDMR: {idmr} does not have a valid origin "
            "(must be paternal or maternal. case insensitive). "
            "This iDMR will not be used for PofO assignment and PofO "
            "assignment score calculation."
        )
    return dmr_chrom, dmr_start, dmr_end, origin


def get_dmr_records(
    tb_dmltest: Any, tb_calldml: Any, dmr_chrom: str, dmr_start: int, dmr_end: int
) -> Tuple[Union[str, List], Union[str, List]]:
    """Query DMR records from tabix files.

    Args:
        tb_dmltest: Tabix iterator for DML test file
        tb_calldml: Tabix iterator for called DML file
        dmr_chrom: Chromosome
        dmr_start: Start position
        dmr_end: End position

    Returns:
        Tuple of (test records, called records)
    """
    try:
        records_all = tb_dmltest.query(dmr_chrom, dmr_start, dmr_end)
    except:
        records_all = "NA"

    try:
        records = tb_calldml.query(dmr_chrom, dmr_start, dmr_end)
    except:
        records = "NA"

    return records_all, records


def process_dmr_methylation(
    records_all: Union[str, List], records: Union[str, List]
) -> Tuple[int, float, float, int, int, float, float]:
    """Process methylation data from DMR records.

    Args:
        records_all: All test records
        records: Called DML records

    Returns:
        Tuple containing:
        - num_cg: Number of CpGs
        - hp1_freq: HP1 methylation frequency
        - hp2_freq: HP2 methylation frequency
        - diff_cg_hp1: Number of HP1 differential CpGs
        - diff_cg_hp2: Number of HP2 differential CpGs
        - diff_cg_hp1_meth: HP1 methylation difference
        - diff_cg_hp2_meth: HP2 methylation difference
    """
    num_cg = 0
    hp1_freq = hp2_freq = 0
    diff_cg_hp1 = diff_cg_hp2 = 0
    diff_cg_hp1_meth = diff_cg_hp2_meth = 0

    if records_all != "NA":
        for record_all in records_all:
            num_cg += 1
            hp1_freq += float(record_all[3])
            hp2_freq += float(record_all[4])

    if records != "NA":
        for record in records:
            diff_cg_meth = float(record[3]) - float(record[4])
            if diff_cg_meth > 0:
                diff_cg_hp1_meth += diff_cg_meth
                diff_cg_hp1 += 1
            elif diff_cg_meth < 0:
                diff_cg_hp2_meth += abs(diff_cg_meth)
                diff_cg_hp2 += 1

    if num_cg > 0:
        hp1_freq = round(hp1_freq / num_cg, 5)
        hp2_freq = round(hp2_freq / num_cg, 5)
    if diff_cg_hp1 > 0:
        diff_cg_hp1_meth = round(diff_cg_hp1_meth / diff_cg_hp1, 5)
    if diff_cg_hp2 > 0:
        diff_cg_hp2_meth = round(diff_cg_hp2_meth / diff_cg_hp2, 5)

    return (
        num_cg,
        hp1_freq,
        hp2_freq,
        diff_cg_hp1,
        diff_cg_hp2,
        diff_cg_hp1_meth,
        diff_cg_hp2_meth,
    )


def write_dmr_output(
    out_meth: TextIO, line: List[str], stats: Tuple, status: str
) -> None:
    """Write DMR statistics to output file.

    Args:
        out_meth: Output file handle
        line: Original DMR line
        stats: Tuple of statistics to write
        status: Status message to append
    """
    # num_cg, diff_cg_hp1, diff_cg_hp2, diff_cg_hp1_meth, diff_cg_hp2_meth, hp1_freq, hp2_freq = stats
    (
        num_cg,
        hp1_freq,
        hp2_freq,
        diff_cg_hp1,
        diff_cg_hp2,
        diff_cg_hp1_meth,
        diff_cg_hp2_meth,
    ) = stats

    tab_line = "\t".join(line)
    out_meth.write(
        f"{tab_line}\t{num_cg}\t{diff_cg_hp1}\t{diff_cg_hp2}\t"
        f"{diff_cg_hp1_meth}\t{diff_cg_hp2_meth}\t{hp1_freq}\t{hp2_freq}\t{status}\n"
    )


def update_origin_counts(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    dmr_chrom: str,
    origin: str,
    num_cg: int,
    diff_scores: Tuple[float, float, int, int],
) -> None:
    """Update origin count statistics.

    Args:
        chrom_hp_origin_count: Dictionary to update
        dmr_chrom: Chromosome
        origin: Parental origin
        num_cg: Number of CpGs
        diff_scores: Tuple of (hp1_score, hp2_score, diff_cg_hp1, diff_cg_hp2)
    """
    diff_cg_hp1_score, diff_cg_hp2_score, diff_cg_hp1, diff_cg_hp2 = diff_scores

    if origin == "maternal":
        update_maternal_counts(chrom_hp_origin_count, dmr_chrom, num_cg, diff_scores)
    elif origin == "paternal":
        update_paternal_counts(chrom_hp_origin_count, dmr_chrom, num_cg, diff_scores)


def PofO_dmr(
    known_dmr: str, out: str, min_cg: int, cpg_difference: float
) -> DefaultDict[Tuple[str, str], DefaultDict[str, int]]:
    """Map differentially methylated CpGs to known iDMRs for PofO assignment.

    Args:
        known_dmr: Path to file containing known imprinted DMRs
        out: Output file prefix
        min_cg: Minimum number of CpGs required in a DMR
        cpg_difference: Minimum required difference in methylation between haplotypes

    Returns:
        Nested defaultdict mapping (chrom, parent_origin) to DMR statistics
    """
    # Initialize files and data structures
    out_meth = init_output_files(out, known_dmr)
    tb_calldml = tabix.open(out + "_Temp_callDML.tsv.gz")
    tb_dmltest = tabix.open(out + "_Temp_DMLtest.tsv.gz")
    chrom_hp_origin_count = defaultdict(lambda: defaultdict(int))

    # Process each DMR
    with openfile(known_dmr) as dmr_file:
        next(dmr_file)  # Skip header
        for line in dmr_file:
            line = line.rstrip().split("\t")
            dmr_chrom, dmr_start, dmr_end, origin = process_dmr_line(line)

            # Get records from tabix files
            records_all, records = get_dmr_records(
                tb_dmltest, tb_calldml, dmr_chrom, dmr_start, dmr_end
            )

            # Process methylation data
            meth_stats = process_dmr_methylation(records_all, records)
            num_cg = meth_stats[0]

            # Check validity conditions and write output
            status = check_dmr_validity(num_cg, meth_stats[3:5], min_cg, cpg_difference)
            if status != "Included":
                write_dmr_output(out_meth, line, meth_stats, status)
                continue

            # Process valid DMR
            write_dmr_output(out_meth, line, meth_stats, "Included")

            # Calculate and update scores
            diff_scores = calculate_diff_scores(meth_stats[3:])
            update_origin_counts(
                chrom_hp_origin_count, dmr_chrom, origin, num_cg, diff_scores
            )

    out_meth.close()
    return chrom_hp_origin_count


def pofo_final_dict(
    chrom_hp_origin_count: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    min_pofo_score: float,
) -> DefaultDict[str, Dict[str, List[Union[str, int, float]]]]:
    """Generate final parent-of-origin assignments from DMR statistics.

    Args:
        chrom_hp_origin_count: Input statistics from PofO_dmr containing counts
            per chromosome/parent
        min_pofo_score: Minimum score threshold for PofO assignment

    Returns:
        Dictionary mapping chromosomes to haplotype origin assignments and statistics
    """
    chrom_hp_origin = defaultdict(dict)
    for key, val in chrom_hp_origin_count.items():
        chrom, origin = key
        hp1_score = val["score_HP1"]
        hp2_score = val["score_HP2"]
        hp1_diff = val["diff_HP1"]
        hp2_diff = val["diff_HP2"]
        hp1_dmr_count = val["dmr_HP1"]
        hp2_dmr_count = val["dmr_HP2"]
        hp1_allcg_count = val["allcg_HP1"]
        hp2_allcg_count = val["allcg_HP2"]
        hp1_alldiffcg_count = val["alldiff_HP1"]
        hp2_alldiffcg_count = val["alldiff_HP2"]
        add_info1 = [
            hp1_score,
            hp2_score,
            hp1_diff,
            hp2_diff,
            hp1_dmr_count,
            hp2_dmr_count,
            hp1_alldiffcg_count,
            hp2_alldiffcg_count,
            hp1_allcg_count,
            hp2_allcg_count,
        ]
        add_info2 = [
            hp2_score,
            hp1_score,
            hp2_diff,
            hp1_diff,
            hp2_dmr_count,
            hp1_dmr_count,
            hp2_alldiffcg_count,
            hp1_alldiffcg_count,
            hp2_allcg_count,
            hp1_allcg_count,
        ]
        if origin.lower() == "maternal":
            if (
                hp1_score > hp2_score
                and hp1_score / (hp1_score + hp2_score) > min_pofo_score
            ):
                chrom_hp_origin[chrom]["HP1"] = ["maternal"] + add_info1
                chrom_hp_origin[chrom]["HP2"] = ["paternal"] + add_info1
            elif (
                hp2_score > hp1_score
                and hp2_score / (hp1_score + hp2_score) > min_pofo_score
            ):
                chrom_hp_origin[chrom]["HP2"] = ["maternal"] + add_info2
                chrom_hp_origin[chrom]["HP1"] = ["paternal"] + add_info2
        elif origin.lower() == "paternal":
            if (
                hp1_score > hp2_score
                and hp1_score / (hp1_score + hp2_score) > min_pofo_score
            ):
                chrom_hp_origin[chrom]["HP1"] = ["paternal"] + add_info1
                chrom_hp_origin[chrom]["HP2"] = ["maternal"] + add_info1
            elif (
                hp2_score > hp1_score
                and hp2_score / (hp1_score + hp2_score) > min_pofo_score
            ):
                chrom_hp_origin[chrom]["HP2"] = ["paternal"] + add_info2
                chrom_hp_origin[chrom]["HP1"] = ["maternal"] + add_info2
    return chrom_hp_origin


def process_cpg_mod_freq(
    input_file: str, output_file: str, is_modkit: bool = False
) -> None:
    """Process CpG modification frequency files from methylation tools.

    Handles output from either modkit or aligned_bam_to_cpg_scores to produce
    standardized CpG modification frequency files.

    Args:
        input_file: Path to input BED file
        output_file: Path to output file
        is_modkit: True if input is from modkit, False if from aligned_bam_to_cpg_scores
    """
    # Write header first
    with open(output_file, "w") as out:
        out.write("Chromosome\tStart\tEnd\tCov\tMod\tFreq\n")

    # Process file content
    with open(input_file) as f:
        with open(output_file, "a") as out:
            for line in f:
                fields = line.strip().split("\t")

                if is_modkit:
                    # For modkit output: check if modification type is 'm'
                    if len(fields) > 3 and fields[3] == "m":
                        chrom, start, end = fields[0:3]
                        cov, mod = fields[9], fields[11]
                        freq = float(mod) / float(cov) if float(cov) != 0 else 0
                        out.write(f"{chrom}\t{start}\t{end}\t{cov}\t{mod}\t{freq}\n")
                else:
                    # For aligned_bam_to_cpg_scores output
                    chrom, start, end = fields[0:3]
                    cov, mod = fields[5], fields[6]
                    freq = float(mod) / float(cov) if float(cov) != 0 else 0
                    out.write(f"{chrom}\t{start}\t{end}\t{cov}\t{mod}\t{freq}\n")


def out_freq_methbam(
    out: str, processes: int, reference: str, pb_tech: bool
) -> None:
    """Run methylation frequency analysis and process outputs.

    Runs either modkit or aligned_bam_to_cpg_scores based on technology,
    then processes the outputs into standardized format.

    Args:
        out: Output file prefix path
        processes: Number of parallel processes to use
        reference: Path to reference genome
        pb_tech: True if using PacBio data, False if ONT
    """
    out_freqhp1 = out + "_Temp_NonPofO_HP1-HP2_MethylationHP1.tsv"
    out_freqhp2 = out + "_Temp_NonPofO_HP1-HP2_MethylationHP2.tsv"
    out_dir = os.path.dirname(out)
    out_pref = os.path.basename(out)
    if pb_tech:

        run_command(
            "aligned_bam_to_cpg_scores --bam {} --output-prefix {}"
            " --threads {} --modsites-mode reference "
            "--ref {}".format(
                out + "_Temp-NonPofO_dmr.bam",
                out + "_Temp-NonPofO_CpGModFreq",
                processes,
                reference,
            ),
        )

        # process cpg for aligned_bam_to_cpg_scores output
        process_cpg_mod_freq(
            out + "_Temp-NonPofO_CpGModFreq.hap1.bed.gz", out_freqhp1, is_modkit=False
        )
        process_cpg_mod_freq(
            out + "_Temp-NonPofO_CpGModFreq.hap2.bed.gz", out_freqhp2, is_modkit=False
        )
    else:
        run_command(
            "modkit pileup -t {} --prefix {} "
            "--partition-tag HP --combine-strands --cpg -r {} "
            "{} {}".format(
                processes,
                out_pref + "_Temp-NonPofO_CpGModFreq",
                reference,
                out + "_Temp-NonPofO_dmr.bam",
                out_dir,
            ),
        )

        # process cpg for modkit output
        process_cpg_mod_freq(
            out + "_Temp-NonPofO_CpGModFreq_1.bed", out_freqhp1, is_modkit=True
        )
        process_cpg_mod_freq(
            out + "_Temp-NonPofO_CpGModFreq_2.bed", out_freqhp2, is_modkit=True
        )

    run_command(
        "rm {}*".format(out + "_Temp-NonPofO_CpGModFreq"),
    )
