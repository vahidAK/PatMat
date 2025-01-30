from collections import defaultdict
from typing import DefaultDict, Dict, Iterator, Set, Tuple


def read_info_iter(read_info_file: str) -> Iterator[Tuple[list, list, Tuple[str, str]]]:
    """Iterate over read information from a variant read info file.

    Args:
        read_info_file: Path to the variant read info file

    Yields:
        Tuple containing:
            - line: List of fields from the input line
            - read_info: List of read information fields
            - read_key: Tuple of (chromosome, read_id)
    """
    with open(read_info_file) as VarReadInfo:
        for line in VarReadInfo:
            line = line.rstrip().split("\t")
            for read_info in line[3].split(","):
                read_info = read_info.split(":")
                yield (line, read_info, (line[0], read_info[0]))


def build_read_dict_HP_temp(
    read_info_file: str,
) -> DefaultDict[Tuple[str, str], DefaultDict[str, int]]:
    """Build dictionary mapping read IDs to their haplotype counts.

    Processes primary reads first, then secondary reads if they weren't already processed.

    Args:
        read_info_file: Path to the variant read info file

    Returns:
        Nested defaultdict mapping (chrom, read_id) -> haplotype -> count
    """
    read_dict_HP_temp = defaultdict(lambda: defaultdict(int))
    added_prim_read = set()
    for line, read_info, read_key in read_info_iter(read_info_file):
        if read_info[-1] != "NA":
            if read_info[1] in ["0", "16"]:
                read_dict_HP_temp[read_key][read_info[-1]] += 1
                added_prim_read.add(read_info[0])
    for line, read_info, read_key in read_info_iter(read_info_file):
        if read_info[-1] != "NA" and read_info[0] not in added_prim_read:
            read_dict_HP_temp[read_key][read_info[-1]] += 1
    return read_dict_HP_temp


def build_variant_dict_HP(
    reads_hap: Dict[Tuple[str, str], int],
    read_info_file: str,
    per_var_info: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
) -> DefaultDict[Tuple[str, str], DefaultDict[Tuple[str, int], int]]:
    """Build dictionary mapping variants to haplotype-specific allele counts.

    Args:
        reads_hap: Dictionary mapping (chrom, read_id) to haplotype assignment
        read_info_file: Path to the variant read info file
        per_var_info: Dictionary for storing per-variant statistics

    Returns:
        Nested defaultdict mapping (chrom, pos) -> (allele, haplotype) -> count
    """
    variant_dict_HP = defaultdict(lambda: defaultdict(int))
    for line, read_info, read_key in read_info_iter(read_info_file):
        var_key = (line[0], line[1])
        hap = reads_hap.get(read_key)
        if hap:
            per_var_info[var_key][f"h{hap}all"] += 1
            if line[2] != "noninfo":
                variant_dict_HP[var_key][(line[2], hap)] += 1
    return variant_dict_HP


def add_reads_hap(
    hapRatio: float,
    minvariant: int,
    reads_hap: Dict[Tuple[str, str], int],
    read_dict_HP_temp_reass: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
) -> None:
    """Add reads to haplotype assignments based on count ratios.

    Args:
        hapRatio: Minimum ratio threshold for haplotype assignment
        minvariant: Minimum variant count required for assignment
        reads_hap: Dictionary to store read haplotype assignments
        read_dict_HP_temp_reass: Dictionary containing reassigned read counts
    """
    for key, val in read_dict_HP_temp_reass.items():
        hp1_count = val["1"]
        hp2_count = val["2"]
        if (
            hp1_count > hp2_count
            and hp1_count / (hp1_count + hp2_count) >= hapRatio
            and hp1_count >= minvariant
        ):
            reads_hap[key] = 1
        elif (
            hp2_count > hp1_count
            and hp2_count / (hp1_count + hp2_count) >= hapRatio
            and hp2_count >= minvariant
        ):
            reads_hap[key] = 2


def update_per_var_info(
    hap_ratio: float,
    min_variant: int,
    read_info_file: str,
    read_dict_HP_temp: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
) -> Tuple[
    DefaultDict[Tuple[str, str], DefaultDict[str, int]],
    DefaultDict[Tuple[str, str], DefaultDict[str, int]],
]:
    """Update per-variant information and reassign read counts.

    Args:
        hap_ratio: Minimum ratio threshold for haplotype assignment
        min_variant: Minimum variant count required for assignment
        read_info_file: Path to the variant read info file
        read_dict_HP_temp: Dictionary containing initial read counts

    Returns:
        Tuple containing:
            - per_var_info: Dictionary mapping variants to their statistics
            - read_dict_HP_temp_reass: Dictionary containing reassigned read counts
    """
    reads_hap_temp = build_reads_hap_temp(hap_ratio, min_variant, read_dict_HP_temp)

    per_var_info = defaultdict(lambda: defaultdict(int))
    read_dict_HP_temp_reass = defaultdict(lambda: defaultdict(int))
    with open(read_info_file) as VarReadInfo:
        for line in VarReadInfo:
            line = line.rstrip().split("\t")
            read_count = hp1_count = hp2_count = hp1_count_read = hp2_count_read = 0
            for read_info in line[3].split(","):
                read_info = read_info.split(":")
                read_key = (line[0], read_info[0])
                per_var_info[(line[0], line[1])]["all"] += 1
                read_count += 1
                hp1_count_read += read_dict_HP_temp[read_key]["1"]
                hp2_count_read += read_dict_HP_temp[read_key]["2"]
                if read_key in reads_hap_temp and line[2] != "noninfo":
                    if reads_hap_temp[read_key] == 1:
                        hp1_count += 1
                    elif reads_hap_temp[read_key] == 2:
                        hp2_count += 1

            per_var_info[tuple(line[0:2])][(line[2], "ave1")] = round(
                hp1_count_read / read_count, 5
            )
            per_var_info[tuple(line[0:2])][(line[2], "ave2")] = round(
                hp2_count_read / read_count, 5
            )
            for read_info in line[3].split(","):
                read_info = read_info.split(":")
                read_key = (line[0], read_info[0])
                if (
                    hp1_count > hp2_count
                    and hp1_count / (hp1_count + hp2_count) >= hap_ratio
                    and hp1_count >= min_variant
                ):
                    read_dict_HP_temp_reass[read_key]["1"] += 1
                elif (
                    hp2_count > hp1_count
                    and hp2_count / (hp1_count + hp2_count) >= hap_ratio
                    and hp2_count >= min_variant
                ):
                    read_dict_HP_temp_reass[read_key]["2"] += 1

    return per_var_info, read_dict_HP_temp_reass


def build_reads_hap_temp(
    hapRatio: float,
    minvariant: int,
    read_dict_HP_temp: DefaultDict[Tuple[str, str], DefaultDict[str, int]],
) -> Dict[Tuple[str, str], int]:
    """Build temporary read haplotype assignments based on count ratios.

    Args:
        hapRatio: Minimum ratio threshold for haplotype assignment
        minvariant: Minimum variant count required for assignment
        read_dict_HP_temp: Dictionary containing read counts

    Returns:
        Dictionary mapping (chrom, read_id) to haplotype assignment (1 or 2)
    """
    reads_hap_temp = dict()
    for key, val in read_dict_HP_temp.items():
        hp1_count = val["1"]
        hp2_count = val["2"]
        if (
            hp1_count > hp2_count
            and hp1_count / (hp1_count + hp2_count) >= hapRatio
            and hp1_count >= minvariant
        ):
            reads_hap_temp[key] = 1
        elif (
            hp2_count > hp1_count
            and hp2_count / (hp1_count + hp2_count) >= hapRatio
            and hp2_count >= minvariant
        ):
            reads_hap_temp[key] = 2
    return reads_hap_temp
