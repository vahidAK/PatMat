def write_merged_dmr_regions(known_dmr_file, output_file):
    """Process DMR regions by extending boundaries and merging overlaps.

    Args:
        known_dmr_file (str): Path to input DMR file
        output_file (str): Path to output processed file
    """
    # Read and process regions
    regions = []
    with open(known_dmr_file) as f:
        next(f)  # Skip header (replaces sed '1d')
        for line in f:
            chrom, start, end = line.split("\t")[:3]
            # Extend regions by 100kb (replaces awk with arithmetic)
            new_start = max(0, int(start) - 100000)  # Ensures start >= 0
            new_end = int(end) + 100000
            regions.append((chrom, new_start, new_end))

    # Sort regions by chromosome and start position
    regions.sort(key=lambda x: (x[0], x[1]))

    # Merge overlapping regions (replaces bedtools merge)
    merged = []
    if not regions:
        return

    current = list(regions[0])

    for region in regions[1:]:
        if (
            region[0] == current[0] and region[1] <= current[2]  # Same chromosome
        ):  # Overlapping regions
            # Extend end position if new end is greater
            current[2] = max(current[2], region[2])
        else:
            merged.append(current)
            current = list(region)
    merged.append(current)

    # Write merged regions
    with open(output_file, "w") as f:
        for region in merged:
            f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")


def write_scores(out, chrom_hp_origin):
    out_scores = open(out + "_PofO_Scores.tsv", "w")
    out_scores.write(
        "Chromosome\tOrigin_HP1\tOrigin_HP2\tPofO_Assignment_Score\t"
        "NormalizedNum_Differentially_Methylated_CGs_Supported_PofO_Assignment\t"
        "NormalizedNum_Differentially_Methylated_CGs_Conflicted_PofO_Assignment\t"
        "Num_Differentially_Methylated_CGs_Supported_PofO_Assignment\t"
        "Num_Differentially_Methylated_CGs_Conflicted_PofO_Assignment\t"
        "Num_iDMRs_Supported_PofO_Assignment\t"
        "Num_iDMRs_Conflicted_PofO_Assignment\t"
        "Num_All_Differentially_Methylated_CGs_At_Supporting_iDMRs\t"
        "Num_All_Differentially_Methylated_CGs_At_Conflicting_iDMRs\t"
        "Num_All_CGs_CouldBeExaminedInBothHaplotypes_At_Supporting_iDMRs\t"
        "Num_All_CGs_CouldBeExaminedInBothHaplotypes_At_Conflicting_iDMRs\n"
    )
    for chrom, val in chrom_hp_origin.items():
        for hp, score in val.items():
            if hp == "HP1":
                origin_hp1 = score[0]
                if origin_hp1 == "maternal":
                    origin_hp2 = "Paternal"
                    origin_hp1 = "Maternal"
                elif origin_hp1 == "paternal":
                    origin_hp2 = "Maternal"
                    origin_hp1 = "Paternal"
                out_scores.write(
                    "\t".join(
                        [
                            chrom,
                            origin_hp1,
                            origin_hp2,
                            str(round(score[1] / (score[1] + score[2]), 5)),
                        ]
                        + list(map(str, score[1:]))
                    )
                    + "\n"
                )
    out_scores.close()
