import os
import subprocess
import warnings
from collections import defaultdict

import tabix

from patmat.io.file_utils import openfile


def PofO_dmr(known_dmr, out, out_meth, min_cg, cpg_difference):
    """
    This function maps differentially methylated CpGs
    to the known iDMRs for PofO assignment
    """
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
            warnings.warn(
                "iDMR {}:{}-{} was ignored because it does not have any CpG in "
                "DMLtest file.".format(dmr_chrom, dmr_start, dmr_end)
            )
            records_all = "NA"
        try:
            records = tb_calldml.query(dmr_chrom, dmr_start, dmr_end)
        except:
            warnings.warn(
                "iDMR {}:{}-{} was ignored because it does not have any CpG in "
                "DMLtest file.".format(dmr_chrom, dmr_start, dmr_end)
            )
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


def out_freq_methbam(out, processes, reference, pbcg, pb_tech):
    out_freqhp1 = out + "_Temp_NonPofO_HP1-HP2_MethylationHP1.tsv"
    out_freqhp2 = out + "_Temp_NonPofO_HP1-HP2_MethylationHP2.tsv"
    out_dir = os.path.dirname(out)
    out_pref = os.path.basename(out)
    if not pb_tech:
        subprocess.run(
            "modkit pileup -t {} --prefix {} "
            "--partition-tag HP --combine-strands --cpg -r {} "
            "{} {}".format(
                processes,
                out_pref + "_Temp-NonPofO_CpGModFreq",
                reference,
                out + "_Temp-NonPofO_dmr.bam",
                out_dir,
            ),
            shell=True,
            check=True,
        )

        subprocess.run(
            "awk -F'\t' '$4==\"m\" {{print $1,$2,$3,$10,$12,$12/$10}}' OFS='\t' {} | "
            "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
            "".format(out + "_Temp-NonPofO_CpGModFreq_1.bed", out_freqhp1),
            shell=True,
            check=True,
        )
        subprocess.run(
            "awk -F'\t' '$4==\"m\" {{print $1,$2,$3,$10,$12,$12/$10}}' OFS='\t' {} | "
            "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
            "".format(out + "_Temp-NonPofO_CpGModFreq_2.bed", out_freqhp2),
            shell=True,
            check=True,
        )
    else:
        subprocess.run(
            "aligned_bam_to_cpg_scores --bam {} --output-prefix {}"
            " --model {} --threads {} --modsites-mode reference "
            "--ref {}".format(
                out + "_Temp-NonPofO_dmr.bam",
                out + "_Temp-NonPofO_CpGModFreq",
                pbcg,
                processes,
                reference,
            ),
            shell=True,
            check=True,
        )

        subprocess.run(
            "awk -F'\t' '{{print $1,$2,$3,$6,$7,$7/$6}}' OFS='\t' {} | "
            "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
            "".format(out + "_Temp-NonPofO_CpGModFreq.hap1.bed", out_freqhp1),
            shell=True,
            check=True,
        )
        subprocess.run(
            "awk -F'\t' '{{print $1,$2,$3,$6,$7,$7/$6}}' OFS='\t' {} | "
            "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
            "".format(out + "_Temp-NonPofO_CpGModFreq.hap2.bed", out_freqhp2),
            shell=True,
            check=True,
        )
    subprocess.run(
        "rm {}*".format(out + "_Temp-NonPofO_CpGModFreq"), shell=True, check=True
    )
