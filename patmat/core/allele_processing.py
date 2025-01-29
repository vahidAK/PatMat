from collections import defaultdict


## ALLELE STATS:
def build_additional_info(
    counts, ratios, line, phase_block_stat=None, blocks_dict=None
):
    """Build additional variant info for output.

    Args:
        counts (dict): Variant allele counts
        ratios (dict): Variant allele ratios
        line (list): VCF line fields
        phase_block_stat (dict, optional): Phase block statistics
        blocks_dict (dict, optional): Block information
    """
    base_info = [
        counts["all_cov"],
        counts["hp1_cov"],
        counts["hp2_cov"],
        counts["hp1_count_ref"],
        counts["hp2_count_ref"],
        counts["hp1_count_alt"],
        counts["hp2_count_alt"],
        ratios["hp1_frac"],
        ratios["hp2_frac"],
        ratios["hp1_frac_ref"],
        ratios["hp2_frac_ref"],
        ratios["hp1_frac_alt"],
        ratios["hp2_frac_alt"],
        ratios["hp1_count_ave_ref"],
        ratios["hp2_count_ave_ref"],
        ratios["hp1_count_ave_alt"],
        ratios["hp2_count_ave_alt"],
    ]

    # Add phase block info if available
    if (
        phase_block_stat
        and blocks_dict
        and (line[0], line[9].split(":")[-1], "agreement") in phase_block_stat
    ):
        block_id = line[9].split(":")[-1]
        block_info = [
            blocks_dict[(line[0], block_id)][0],
            blocks_dict[(line[0], block_id)][1],
            phase_block_stat[(line[0], block_id, "agreement")],
            phase_block_stat[(line[0], block_id, "disagreement")],
        ]
    else:
        block_info = ["NA", "NA", "NA", "NA"]

    return list(map(str, base_info + block_info))


def determine_reassignment(
    counts,
    ratios,
    format_values,
    min_read_reassignment,
    reassignment_format,
):
    """Determine variant reassignment based on allele stats.

    Args:
        counts (dict): Variant allele counts
        ratios (dict): Variant allele ratios
        format_values (dict): Updated format fields
        min_read_reassignment (int): Minimum reads for reassignment
        reassignment_format (tuple): Format strings for reassignment (hp1, hp2, label1, label2)
    """
    hp1_fmt, hp2_fmt, label1, label2 = reassignment_format

    # Check conditions for HP2|HP1 reassignment
    if (
        counts["hp1_count_alt"] > counts["hp2_count_alt"]
        and ratios["hp1_alt_ratio"] > ratios["hp2_alt_ratio"]
        and ratios["hp1_alt_ratio"] >= ratios["hp1_ref_ratio"]
        and counts["hp1_count_alt"] >= min_read_reassignment
    ) or (
        counts["hp2_count_ref"] > counts["hp1_count_ref"]
        and ratios["hp2_ref_ratio"] > ratios["hp1_ref_ratio"]
        and ratios["hp2_ref_ratio"] >= ratios["hp2_alt_ratio"]
        and counts["hp2_count_ref"] >= min_read_reassignment
    ):
        format_values["GT"] = hp1_fmt
        format_values["PS"] = label1

    # Check conditions for HP1|HP2 reassignment
    elif (
        counts["hp2_count_alt"] > counts["hp1_count_alt"]
        and ratios["hp2_alt_ratio"] > ratios["hp1_alt_ratio"]
        and ratios["hp2_alt_ratio"] >= ratios["hp2_ref_ratio"]
        and counts["hp2_count_alt"] >= min_read_reassignment
    ) or (
        counts["hp1_count_ref"] > counts["hp2_count_ref"]
        and ratios["hp1_ref_ratio"] > ratios["hp2_ref_ratio"]
        and ratios["hp1_ref_ratio"] >= ratios["hp1_alt_ratio"]
        and counts["hp1_count_ref"] >= min_read_reassignment
    ):
        format_values["GT"] = hp2_fmt
        format_values["PS"] = label2

    # Default - no reassignment
    else:
        format_values["GT"] = (
            format_values["GT"]
            .replace("|", "/")
            .replace("2/1", "1/2")
            .replace("1/0", "0/1")
        )


def calculate_allele_stats(
    variant_dict_HP, per_var_info, var_key, ref_allele, alt_allele
):
    """Calculate allele counts and ratios for a variant."""
    counts = {
        "hp1_count_alt": variant_dict_HP[var_key][(alt_allele, 1)],
        "hp2_count_alt": variant_dict_HP[var_key][(alt_allele, 2)],
        "hp1_count_ref": variant_dict_HP[var_key][(ref_allele, 1)],
        "hp2_count_ref": variant_dict_HP[var_key][(ref_allele, 2)],
        "hp1_cov": per_var_info[var_key]["h1all"],
        "hp2_cov": per_var_info[var_key]["h2all"],
        "all_cov": per_var_info[var_key]["all"],
    }

    ratios = calculate_coverage_ratios(
        counts, per_var_info, var_key, ref_allele, alt_allele
    )
    return counts, ratios


def calculate_coverage_ratios(counts, per_var_info, var_key, ref_allele, alt_allele):
    """Calculate coverage ratios for alleles."""
    ratios = defaultdict(int)  # Default all ratios to 0

    if counts["all_cov"] > 0:
        ratios["hp1_frac"] = round(counts["hp1_cov"] / counts["all_cov"], 5)
        ratios["hp2_frac"] = round(counts["hp2_cov"] / counts["all_cov"], 5)

    if counts["hp1_cov"] > 0:
        ratios.update(
            {
                "hp1_count_ave_ref": per_var_info[var_key][(ref_allele, "ave1")],
                "hp1_count_ave_alt": per_var_info[var_key][(alt_allele, "ave1")],
                "hp1_frac_ref": round(counts["hp1_count_ref"] / counts["hp1_cov"], 5),
                "hp1_frac_alt": round(counts["hp1_count_alt"] / counts["hp1_cov"], 5),
            }
        )
        total = counts["hp1_count_alt"] + counts["hp1_count_ref"]
        if total > 0:
            ratios["hp1_alt_ratio"] = counts["hp1_count_alt"] / total
            ratios["hp1_ref_ratio"] = counts["hp1_count_ref"] / total

    if counts["hp2_cov"] > 0:
        ratios.update(
            {
                "hp2_count_ave_ref": per_var_info[var_key][(ref_allele, "ave2")],
                "hp2_count_ave_alt": per_var_info[var_key][(alt_allele, "ave2")],
                "hp2_frac_ref": round(counts["hp2_count_ref"] / counts["hp2_cov"], 5),
                "hp2_frac_alt": round(counts["hp2_count_alt"] / counts["hp2_cov"], 5),
            }
        )
        total = counts["hp2_count_alt"] + counts["hp2_count_ref"]
        if total > 0:
            ratios["hp2_alt_ratio"] = counts["hp2_count_alt"] / total
            ratios["hp2_ref_ratio"] = counts["hp2_count_ref"] / total

    return ratios


def process_variant(
    line,
    format_values,
    variant_dict_HP,
    per_var_info,
    min_read_reassignment,
    phase_block_stat=None,
    blocks_dict=None,
):
    """Process a single variant, either biallelic or multiallelic."""
    var_key = tuple(line[0:2])

    # Handle biallelic vs multiallelic
    if format_values["GT"] in ("0/1", "1/0", "0|1", "1|0"):
        ref_allele = line[3].upper()
        alt_allele = line[4].upper()
        reassignment_format = ("1|0", "0|1", "HP2|HP1", "HP1|HP2")
    else:  # Multiallelic
        alt_alleles = line[4].split(",")
        ref_allele = alt_alleles[0].upper()
        alt_allele = alt_alleles[1].upper()
        reassignment_format = ("1|2", "1|2", "Ref_HP2|HP1", "Ref_HP1|HP2")

    counts, ratios = calculate_allele_stats(
        variant_dict_HP, per_var_info, var_key, ref_allele, alt_allele
    )
    additional_info = build_additional_info(
        counts, ratios, line, phase_block_stat, blocks_dict
    )

    determine_reassignment(
        counts,
        ratios,
        format_values,
        min_read_reassignment,
        reassignment_format,
    )

    # Return variant columns
    return (
        line[0:8]
        + [":".join(format_values.keys())]
        + [":".join(format_values.values())]
        + additional_info
    )
