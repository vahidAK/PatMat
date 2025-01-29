import os
import subprocess
import warnings
from collections import defaultdict

import tabix

from patmat.core.subprocess_calls import run_command
from patmat.io.file_utils import openfile


def PofO_dmr(known_dmr, out, min_cg, cpg_difference):
    """
    This function maps differentially methylated CpGs
    to the known iDMRs for PofO assignment
    """

    out_meth = open(out + "_CpG-Methylation-Status-at-DMRs.tsv", "w")
    dmr_file = openfile(known_dmr)
    header = next(dmr_file).rstrip()
    out_meth.write(
        header + "\tAll_CpGs_At_iDMR_CouldBeExaminedInBothHaplotypes\t"
        "DifferentiallyMethylatedCpGs_HypermethylatedOnHP1\t"
        "DifferentiallyMethylatedCpGs_HypermethylatedOnHP2\t"
        "MeanDiffMethylationOf_DifferentiallyMethylatedCpGs_HypermethylatedOnHP1\t"
        "MeanDiffMethylationOf_DifferentiallyMethylatedCpGs_HypermethylatedOnHP2\t"
        "MethylationFrequency_HP1\tMethylationFrequency_HP2\t"
        "Included_Or_Ignored_For_PofO_Assignment\n"
    )
    dmr_file.close()

    tb_calldml = tabix.open(out + "_Temp_callDML.tsv.gz")
    tb_dmltest = tabix.open(out + "_Temp_DMLtest.tsv.gz")
    dmr_file = openfile(known_dmr)
    header = next(dmr_file).rstrip()
    chrom_hp_origin_count = defaultdict(lambda: defaultdict(int))
    for line in dmr_file:
        line = line.rstrip().split("\t")
        dmr_chrom = line[0]
        dmr_start = int(line[1]) - 1
        dmr_end = int(line[2])
        origin = line[3].lower()
        if origin not in ["maternal", "paternal"]:
            warnings.warn(
                "iDMR: {} does not have a valid origin "
                "(must be paternal or maternal. case insensitive). "
                "This iDMR will not be used for PofO assignment and PofO"
                " assignment score calculation.".format("\t".join(line))
            )
        try:
            records_all = tb_dmltest.query(dmr_chrom, dmr_start, dmr_end)
        except:
            # warnings.warn(
            #     "iDMR {}:{}-{} was ignored because it does not have any CpG in "
            #     "DMLtest file.".format(dmr_chrom, dmr_start, dmr_end)
            # )
            records_all = "NA"
        try:
            records = tb_calldml.query(dmr_chrom, dmr_start, dmr_end)
        except:
            # warnings.warn(
            #     "iDMR {}:{}-{} was ignored because it does not have any CpG in "
            #     "DMLtest file.".format(dmr_chrom, dmr_start, dmr_end)
            # )
            records = "NA"
        num_cg = 0
        hp1_freq = 0
        hp2_freq = 0
        diff_cg_hp1 = 0
        diff_cg_hp2 = 0
        diff_cg_hp1_meth = 0
        diff_cg_hp2_meth = 0
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
        if num_cg < 1:
            out_meth.write(
                "\t".join(line)
                + "\t"
                + str(num_cg)
                + "\tNA\tNA\tNA\tNA\tNA\tNA\tIgnored:No_CpG\n"
            )
            continue
        hp1_freq = round(hp1_freq / num_cg, 5)
        hp2_freq = round(hp2_freq / num_cg, 5)
        if diff_cg_hp1 > 0:
            diff_cg_hp1_meth = round(diff_cg_hp1_meth / diff_cg_hp1, 5)
        else:
            diff_cg_hp1_meth = 0
        if diff_cg_hp2 > 0:
            diff_cg_hp2_meth = round(diff_cg_hp2_meth / diff_cg_hp2, 5)
        else:
            diff_cg_hp2_meth = 0
        if num_cg < min_cg and abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
            out_meth.write(
                "\t".join(line)
                + "\t"
                + str(num_cg)
                + "\t"
                + str(diff_cg_hp1)
                + "\t"
                + str(diff_cg_hp2)
                + "\t"
                + str(diff_cg_hp1_meth)
                + "\t"
                + str(diff_cg_hp2_meth)
                + "\t"
                + str(hp1_freq)
                + "\t"
                + str(hp2_freq)
                + "\t"
                + "Ignored:DidNotMeet_min_cg_And_cpg_difference\n"
            )
            continue
        elif num_cg < min_cg and diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
            out_meth.write(
                "\t".join(line)
                + "\t"
                + str(num_cg)
                + "\t"
                + str(diff_cg_hp1)
                + "\t"
                + str(diff_cg_hp2)
                + "\t"
                + str(diff_cg_hp1_meth)
                + "\t"
                + str(diff_cg_hp2_meth)
                + "\t"
                + str(hp1_freq)
                + "\t"
                + str(hp2_freq)
                + "\t"
                + "Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylated"
                "CpGsIsZeroInBothHaplotypes\n"
            )
            continue
        elif num_cg < min_cg:
            out_meth.write(
                "\t".join(line)
                + "\t"
                + str(num_cg)
                + "\t"
                + str(diff_cg_hp1)
                + "\t"
                + str(diff_cg_hp2)
                + "\t"
                + str(diff_cg_hp1_meth)
                + "\t"
                + str(diff_cg_hp2_meth)
                + "\t"
                + str(hp1_freq)
                + "\t"
                + str(hp2_freq)
                + "\t"
                + "Ignored:DidNotMeet_min_cg\n"
            )
            continue
        elif abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
            out_meth.write(
                "\t".join(line)
                + "\t"
                + str(num_cg)
                + "\t"
                + str(diff_cg_hp1)
                + "\t"
                + str(diff_cg_hp2)
                + "\t"
                + str(diff_cg_hp1_meth)
                + "\t"
                + str(diff_cg_hp2_meth)
                + "\t"
                + str(hp1_freq)
                + "\t"
                + str(hp2_freq)
                + "\t"
                + "Ignored:DidNotMeet_cpg_difference\n"
            )
            continue
        elif diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
            out_meth.write(
                "\t".join(line)
                + "\t"
                + str(num_cg)
                + "\t"
                + str(diff_cg_hp1)
                + "\t"
                + str(diff_cg_hp2)
                + "\t"
                + str(diff_cg_hp1_meth)
                + "\t"
                + str(diff_cg_hp2_meth)
                + "\t"
                + str(hp1_freq)
                + "\t"
                + str(hp2_freq)
                + "\t"
                + "Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes\n"
            )
            continue
        out_meth.write(
            "\t".join(line)
            + "\t"
            + str(num_cg)
            + "\t"
            + str(diff_cg_hp1)
            + "\t"
            + str(diff_cg_hp2)
            + "\t"
            + str(diff_cg_hp1_meth)
            + "\t"
            + str(diff_cg_hp2_meth)
            + "\t"
            + str(hp1_freq)
            + "\t"
            + str(hp2_freq)
            + "\t"
            + "Included\n"
        )
        diff_cg_hp1_score = diff_cg_hp1 * diff_cg_hp1_meth
        diff_cg_hp2_score = diff_cg_hp2 * diff_cg_hp2_meth
        if origin == "maternal":
            chrom_hp_origin_count[(dmr_chrom, "maternal")][
                "score_HP1"
            ] += diff_cg_hp1_score
            chrom_hp_origin_count[(dmr_chrom, "maternal")][
                "score_HP2"
            ] += diff_cg_hp2_score
            chrom_hp_origin_count[(dmr_chrom, "paternal")][
                "score_HP2"
            ] += diff_cg_hp1_score
            chrom_hp_origin_count[(dmr_chrom, "paternal")][
                "score_HP1"
            ] += diff_cg_hp2_score
            chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP1"] += diff_cg_hp1
            chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP2"] += diff_cg_hp2
            chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP2"] += diff_cg_hp1
            chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP1"] += diff_cg_hp2
            if diff_cg_hp1_score > diff_cg_hp2_score:
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
            elif diff_cg_hp1_score < diff_cg_hp2_score:
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

        elif origin == "paternal":
            chrom_hp_origin_count[(dmr_chrom, "maternal")][
                "score_HP1"
            ] += diff_cg_hp2_score
            chrom_hp_origin_count[(dmr_chrom, "maternal")][
                "score_HP2"
            ] += diff_cg_hp1_score
            chrom_hp_origin_count[(dmr_chrom, "paternal")][
                "score_HP2"
            ] += diff_cg_hp2_score
            chrom_hp_origin_count[(dmr_chrom, "paternal")][
                "score_HP1"
            ] += diff_cg_hp1_score
            chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP1"] += diff_cg_hp2
            chrom_hp_origin_count[(dmr_chrom, "maternal")]["diff_HP2"] += diff_cg_hp1
            chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP2"] += diff_cg_hp2
            chrom_hp_origin_count[(dmr_chrom, "paternal")]["diff_HP1"] += diff_cg_hp1
            if diff_cg_hp1_score > diff_cg_hp2_score:
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
            elif diff_cg_hp1_score < diff_cg_hp2_score:
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

    dmr_file.close()
    out_meth.close()
    return chrom_hp_origin_count


def pofo_final_dict(chrom_hp_origin_count, min_pofo_score):
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


def process_cpg_mod_freq(input_file, output_file, is_modkit=False):
    """Process CpG modification frequency files from either modkit or aligned_bam_to_cpg_scores.

    Args:
        input_file (str): Path to input bed file
        output_file (str): Path to output file
        is_modkit (bool): True if input is from modkit, False if from aligned_bam_to_cpg_scores
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


def out_freq_methbam(out, processes, reference, pbcg, pb_tech):
    out_freqhp1 = out + "_Temp_NonPofO_HP1-HP2_MethylationHP1.tsv"
    out_freqhp2 = out + "_Temp_NonPofO_HP1-HP2_MethylationHP2.tsv"
    out_dir = os.path.dirname(out)
    out_pref = os.path.basename(out)
    if pb_tech:

        run_command(
            "aligned_bam_to_cpg_scores --bam {} --output-prefix {}"
            " --model {} --threads {} --modsites-mode reference "
            "--ref {}".format(
                out + "_Temp-NonPofO_dmr.bam",
                out + "_Temp-NonPofO_CpGModFreq",
                pbcg,
                processes,
                reference,
            ),
        )

        # process cpg for aligned_bam_to_cpg_scores output
        process_cpg_mod_freq(
            out + "_Temp-NonPofO_CpGModFreq.hap1.bed", out_freqhp1, is_modkit=False
        )
        process_cpg_mod_freq(
            out + "_Temp-NonPofO_CpGModFreq.hap2.bed", out_freqhp2, is_modkit=False
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
