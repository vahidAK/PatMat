from pathlib import Path
from typing import Dict, List, Tuple, Union


def write_merged_dmr_regions(
    known_dmr_file: Union[str, Path],
    output_file: Union[str, Path],
    extension_size: int = 100000,
) -> None:
    """Process DMR regions by extending boundaries and merging overlaps.

    Args:
        known_dmr_file: Path to input DMR file containing regions in BED format (chr, start, end)
        output_file: Path where processed regions will be written
        extension_size: Size in bp to extend regions on each side (default: 100kb)

    Notes:
        - Input file should be tab-separated with header
        - Uses first 3 columns (chromosome, start, end)
        - Extends regions by extension_size on both sides
        - Merges overlapping regions after extension
        - Ensures start positions are not negative
        - Outputs sorted, merged regions in BED format

    Example output format:
        chr1    1000    2000
        chr1    2500    3500
        chr2    1000    2000
    """
    # Read and process regions
    regions: List[Tuple[str, int, int]] = []
    with open(known_dmr_file) as f:
        next(f)  # Skip header
        for line in f:
            chrom, start, end = line.split("\t")[:3]
            new_start = max(0, int(start) - extension_size)
            new_end = int(end) + extension_size
            regions.append((chrom, new_start, new_end))

    if not regions:
        return

    # Sort regions by chromosome and start position
    regions.sort(key=lambda x: (x[0], x[1]))

    # Merge overlapping regions
    merged: List[List[Union[str, int]]] = []
    current = list(regions[0])

    for region in regions[1:]:
        if region[0] == current[0] and region[1] <= current[2]:
            current[2] = max(current[2], region[2])
        else:
            merged.append(current)
            current = list(region)
    merged.append(current)

    # Write merged regions
    with open(output_file, "w") as f:
        for region in merged:
            f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")


def write_scores(
    out: Union[str, Path],
    chrom_hp_origin: Dict[str, Dict[str, List[Union[str, int, float]]]],
) -> None:
    """Write parent-of-origin scores to file.

    Args:
        out: Output file prefix (will append '_PofO_Scores.tsv')
        chrom_hp_origin: Dictionary containing PofO assignments and statistics.
            Structure: {chrom: {'HP1': [origin, stats...], 'HP2': [origin, stats...]}}

    Notes:
        - Only processes HP1 entries for each chromosome
        - Capitalizes maternal/paternal in output
        - Calculates PofO assignment score as score1/(score1 + score2)
        - Rounds score to 5 decimal places

    Output columns:
        - Chromosome
        - Origin_HP1 (Maternal/Paternal)
        - Origin_HP2 (Maternal/Paternal)
        - PofO_Assignment_Score
        - Various methylation and DMR statistics
    """
    with open(f"{out}_PofO_Scores.tsv", "w") as out_scores:
        # Write header
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

        # Process each chromosome
        for chrom, val in chrom_hp_origin.items():
            for hp, score in val.items():
                if hp != "HP1":
                    continue

                # Format origin strings
                origin_hp1, origin_hp2 = format_origins(score[0])

                # Calculate and write scores
                pofo_score = round(score[1] / (score[1] + score[2]), 5)
                scores = "\t".join(map(str, score[1:]))
                out_scores.write(
                    f"{chrom}\t{origin_hp1}\t{origin_hp2}\t{pofo_score}\t" f"{scores}\n"
                )


def format_origins(origin: str) -> Tuple[str, str]:
    """Format origin strings with proper capitalization."""
    if origin == "maternal":
        return "Maternal", "Paternal"
    return "Paternal", "Maternal"
