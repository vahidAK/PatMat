from collections import defaultdict
from typing import DefaultDict, Dict, Optional, Sequence, TextIO, Tuple, Union

from patmat.io.file_utils import openfile


##### Process variant assignments
def process_variant_assignments(
    vcf_file: str,
    out_prefix: str,
    re_assignment_vars: Dict[Tuple[str, str], Sequence[Union[str, int, float]]],
    chrom_hp_origin: Dict[str, Dict[str, Sequence[Union[str, int, float]]]],
) -> DefaultDict[str, DefaultDict[str, int]]:
    """Process variant assignments and write output files.

    Args:
        vcf_file: Path to input VCF file
        out_prefix: Prefix for output files
        re_assignment_vars: Dictionary mapping (chrom, pos) to variant info
        chrom_hp_origin: Dictionary mapping chromosomes to parent-of-origin assignments

    Returns:
        Nested defaultdict tracking variant counts per chromosome
    """

    info_out_dict = defaultdict(lambda: defaultdict(int))

    with openfile(vcf_file) as vf:
        with open(f"{out_prefix}_PofO_Assigned.vcf", "w") as assignment_file:
            with open(f"{out_prefix}_Variant_Assignment_info.tsv", "w") as info_file:
                # Write headers and get first data line
                first_data_line = write_headers(vf, assignment_file, info_file)

                # Process first line if exists
                if first_data_line:
                    line = first_data_line.rstrip().split("\t")
                    process_variant_line(
                        line,
                        re_assignment_vars,
                        chrom_hp_origin,
                        assignment_file,
                        info_file,
                        info_out_dict,
                    )

                # Process remaining variants
                for line in vf:
                    if line.startswith("#"):
                        continue

                    line = line.rstrip().split("\t")
                    process_variant_line(
                        line,
                        re_assignment_vars,
                        chrom_hp_origin,
                        assignment_file,
                        info_file,
                        info_out_dict,
                    )
    return info_out_dict


def write_headers(
    vcf_file: TextIO, assignment_file: TextIO, info_file: TextIO
) -> Optional[str]:
    """Write headers to output files and return first data line if found.

    Args:
        vcf_file: Input VCF file handle
        assignment_file: Main output file handle
        info_file: Info output file handle

    Returns:
        First non-header line if found, None otherwise
    """
    first_data_line = None
    for line in vcf_file:
        if not line.startswith("#"):
            first_data_line = line
            break

        assignment_file.write(line)
        if line.startswith("#CHROM"):
            info_file.write(line.rstrip() + VARIANT_INFO_HEADER)

    return first_data_line


def process_variant_line(
    line: Sequence[str],
    re_assignment_vars: Dict[Tuple[str, str], Sequence[Union[str, int, float]]],
    chrom_hp_origin: Dict[str, Dict[str, Sequence[Union[str, int, float]]]],
    assignment_file: TextIO,
    info_file: TextIO,
    info_out_dict: DefaultDict[str, DefaultDict[str, int]],
) -> None:
    """Process a single variant line.

    Args:
        line: List of VCF fields
        re_assignment_vars: Dictionary mapping (chrom, pos) to variant info
        chrom_hp_origin: Dictionary mapping chromosomes to parent-of-origin assignments
        assignment_file: Main output file handle
        info_file: Info output file handle
        info_out_dict: Dictionary tracking variant counts
    """
    # Check variant type
    var_type = determine_variant_type(line)
    if var_type:
        update_variant_counts(info_out_dict, line[0], var_type)

    # Process reassignment
    if tuple(line[0:2]) in re_assignment_vars:
        process_reassigned_variant(
            line,
            re_assignment_vars,
            chrom_hp_origin,
            var_type,
            assignment_file,
            info_file,
            info_out_dict,
        )
    else:
        process_unassigned_variant(line, assignment_file, info_file)


def determine_variant_type(line: Sequence[str]) -> Optional[str]:
    """Determine if variant is SNV, indel, or neither.

    Args:
        line: List of VCF fields

    Returns:
        'snv', 'indel', or None if not a variant
    """
    if not line[9].startswith(("0/1", "1/0", "0|1", "1|0", "1/2", "1|2", "2/1", "2|1")):
        return None

    if (len(line[3]) == 1 and len(line[4]) == 1) or (
        len(line[3]) == 1 and len(line[4]) == 3 and "," in line[4]
    ):
        return "snv"
    return "indel"


def update_variant_counts(
    info_dict: DefaultDict[str, DefaultDict[str, int]], chrom: str, var_type: str
) -> None:
    """Update variant counts in info dictionary.

    Args:
        info_dict: Dictionary tracking variant counts
        chrom: Chromosome name
        var_type: Type of variant ('snv' or 'indel')
    """
    if var_type == "snv":
        info_dict[chrom]["all_het_snvs"] += 1
    else:
        info_dict[chrom]["all_het_indels"] += 1


def process_reassigned_variant(
    line: Sequence[str],
    re_assignment_vars: Dict[Tuple[str, str], Sequence[Union[str, int, float]]],
    chrom_hp_origin: Dict[str, Dict[str, Sequence[Union[str, int, float]]]],
    var_type: Optional[str],
    assignment_file: TextIO,
    info_file: TextIO,
    info_out_dict: DefaultDict[str, DefaultDict[str, int]],
) -> None:
    """Process a variant that has been reassigned.

    Args:
        line: List of VCF fields
        re_assignment_vars: Dictionary mapping (chrom, pos) to variant info
        chrom_hp_origin: Dictionary mapping chromosomes to parent-of-origin assignments
        var_type: Type of variant ('snv', 'indel', or None)
        assignment_file: Main output file handle
        info_file: Info output file handle
        info_out_dict: Dictionary tracking variant counts
    """

    var_info = re_assignment_vars[tuple(line[0:2])]

    if not is_valid_reassignment(var_info, line, chrom_hp_origin):
        write_unphased_variant(var_info, assignment_file, info_file)
        return

    # Extract counts and update statistics
    counts = extract_variant_counts(var_info)
    update_pofo_counts(info_out_dict, line[0], var_type, var_info)

    # Write output based on parent of origin
    if chrom_hp_origin[line[0]]["HP1"][0] == "maternal":
        write_maternal_variant(var_info, counts, assignment_file, info_file)
    elif chrom_hp_origin[line[0]]["HP1"][0] == "paternal":
        write_paternal_variant(var_info, counts, assignment_file, info_file)


def process_unassigned_variant(
    line: Sequence[str], assignment_file: TextIO, info_file: TextIO
) -> None:
    """Process a variant that hasn't been reassigned.

    Args:
        line: List of VCF fields
        assignment_file: Main output file handle
        info_file: Info output file handle
    """
    format_values = process_format_fields(line)
    out_line = format_unassigned_variant(line, format_values)

    assignment_file.write(out_line + "\n")
    info_file.write(out_line + "\t" + "\t".join(["NA"] * 33) + "\n")


# Constants
VARIANT_INFO_HEADER = (
    "\tNumAllReads\tNumReadsHP1\tNumReadsHP2\t"
    "NumReadsHP1RefOrLeftAllele\tNumReadsHP2RefOrLeftAllele\t"
    "NumReadsHP1AltOrRightAllele\tNumReadsHP2AltOrRightAllele\t"
    "FracReadsHP1\tFracReadsHP2\t"
    "FracHP1ReadsWithRefOrLeftAllele\tFracHP2ReadsWithRefOrLeftAllele\t"
    "FracHP1ReadsWithAltOrRightAllele\tFracHP2ReadsWithAltOrLeftAllele\t"
    "MeanNumherOfHP1-InitialPhasedVariantsAccrossReadsMappedToRef/LeftAllele\t"
    "MeanNumherOfHP2-InitialPhasedVariantsAccrossReadsMappedToRef/LeftAllele\t"
    "MeanNumherOfHP1-InitialPhasedVariantsAccrossReadsMappedToAlt/RightAllele\t"
    "MeanNumherOfHP2-InitialPhasedVariantsAccrossReadsMappedToAlt/RightAllele\t"
    "BlockStart\tBlockEnd\t"
    "NumberOfSupportiveStrandSeqPhasedVariantsAtThePhasedBlock\t"
    "NumberOfConflictingStrandSeqPhasedVariantsAtThePhasedBlock\t"
    "NumReadsMaternal\tNumReadsPaternal\t"
    "NumReadsMaternalRefOrLeftAllele\t"
    "NumReads_PaternalRefOrLeftAllele\t"
    "NumReads_MaternalAltOrRightAllele\t"
    "NumReads_PaternalAltOrRightAllele\t\t"
    "FracReadsMaternal\tFracReadsPaternal\t"
    "FracMaternalReadsWithRefOrLeftAllele\tFracPaternalReadsWithRefOrLeftAllele\t"
    "FracMaternalReadsWithAltOrRightAllele\tFracPaternalReadsWithAltOrRightAllele\n"
)


def is_valid_reassignment(
    var_info: Sequence[Union[str, int, float]],
    line: Sequence[str],
    chrom_hp_origin: Dict[str, Dict[str, Sequence[Union[str, int, float]]]],
) -> bool:
    """Check if variant reassignment is valid.

    Args:
        var_info: List of variant information
        line: List of VCF fields
        chrom_hp_origin: Dictionary mapping chromosomes to parent-of-origin assignments

    Returns:
        True if reassignment is valid, False otherwise
    """
    return var_info[9].startswith(("1|0", "0|1", "1|2")) and line[0] in chrom_hp_origin


def extract_variant_counts(
    var_info: Sequence[Union[str, int, float]]
) -> Dict[str, Union[str, int, float]]:
    """Extract count information from variant info.

    Args:
        var_info: List of variant information

    Returns:
        Dictionary mapping count types to their values
    """
    return {
        "hp1_count": var_info[11],
        "hp2_count": var_info[12],
        "hp1_count_ref": var_info[13],
        "hp2_count_ref": var_info[14],
        "hp1_count_alt": var_info[15],
        "hp2_count_alt": var_info[16],
        "hp1_frac": var_info[17],
        "hp2_frac": var_info[18],
        "hp1_ref_frac": var_info[19],
        "hp2_ref_frac": var_info[20],
        "hp1_alt_frac": var_info[21],
        "hp2_alt_frac": var_info[22],
    }


def update_pofo_counts(
    info_dict: DefaultDict[str, DefaultDict[str, int]],
    chrom: str,
    var_type: Optional[str],
    var_info: Sequence[Union[str, int, float]],
) -> None:
    """Update parent-of-origin counts.

    Args:
        info_dict: Dictionary tracking variant counts
        chrom: Chromosome name
        var_type: Type of variant ('snv', 'indel', or None)
        var_info: List of variant information
    """
    if var_type == "snv" and var_info[9].startswith(("1|0", "0|1", "1|2")):
        info_dict[chrom]["pofo_het_snvs"] += 1
    elif var_type == "indel" and var_info[9].startswith(("1|0", "0|1", "1|2")):
        info_dict[chrom]["pofo_het_indels"] += 1


def write_unphased_variant(
    var_info: Sequence[Union[str, int, float]],
    assignment_file: TextIO,
    info_file: TextIO,
) -> None:
    """Write an unphased variant to output files.

    Args:
        var_info: List of variant information fields
        assignment_file: Main output file handle
        info_file: Info output file handle

    Note:
        Replaces phased separator '|' with unphased '/' and adjusts allele order
    """
    out_line = (
        "\t".join(var_info[0:10])
        .replace(":PS", "")
        .replace("1|0", "0/1")
        .replace("2|1", "1/2")
        .replace("|", "/")
        .replace(":Ref_HP1/HP2", "")
        .replace(":Ref_HP2/HP1", "")
        .replace(":HP1/HP2", "")
        .replace(":HP2/HP1", "")
    )
    assignment_file.write(out_line + "\n")
    info_file.write("\t".join(var_info + ["NA"] * 12) + "\n")


def write_maternal_variant(
    var_info: Sequence[Union[str, int, float]],
    counts: Dict[str, Union[str, int, float]],
    assignment_file: TextIO,
    info_file: TextIO,
) -> None:
    """Write a maternal variant to output files.

    Args:
        var_info: List of variant information fields
        counts: Dictionary of count statistics
        assignment_file: Main output file handle
        info_file: Info output file handle
    """
    out_line = "\t".join(var_info[0:10]).replace("HP1", "Mat").replace("HP2", "Pat")
    assignment_file.write(out_line + "\n")

    count_info = [
        counts["hp1_count"],
        counts["hp2_count"],
        counts["hp1_count_ref"],
        counts["hp2_count_ref"],
        counts["hp1_count_alt"],
        counts["hp2_count_alt"],
        counts["hp1_frac"],
        counts["hp2_frac"],
        counts["hp1_ref_frac"],
        counts["hp2_ref_frac"],
        counts["hp1_alt_frac"],
        counts["hp2_alt_frac"],
    ]

    info_file.write(out_line + "\t" + "\t".join(var_info[10:] + count_info) + "\n")


def write_paternal_variant(
    var_info: Sequence[Union[str, int, float]],
    counts: Dict[str, Union[str, int, float]],
    assignment_file: TextIO,
    info_file: TextIO,
) -> None:
    """Write a paternal variant to output files.

    Args:
        var_info: List of variant information fields
        counts: Dictionary of count statistics
        assignment_file: Main output file handle
        info_file: Info output file handle

    Note:
        Handles different format cases (1|0, 0|1) and adjusts labels accordingly
    """
    if var_info[9].startswith("1|0"):
        out_line = (
            "\t".join(var_info[0:10])
            .replace("1|0", "0|1")
            .replace("HP1", "Pat")
            .replace("HP2", "Mat")
        )
    elif var_info[9].startswith("0|1"):
        out_line = (
            "\t".join(var_info[0:10])
            .replace("0|1", "1|0")
            .replace("HP1", "Pat")
            .replace("HP2", "Mat")
        )
    else:
        out_line = "\t".join(var_info[0:10]).replace("HP1", "Pat").replace("HP2", "Mat")

    assignment_file.write(out_line + "\n")

    count_info = [
        counts["hp2_count"],
        counts["hp1_count"],
        counts["hp2_count_ref"],
        counts["hp1_count_ref"],
        counts["hp2_count_alt"],
        counts["hp1_count_alt"],
        counts["hp2_frac"],
        counts["hp1_frac"],
        counts["hp2_ref_frac"],
        counts["hp1_ref_frac"],
        counts["hp2_alt_frac"],
        counts["hp1_alt_frac"],
    ]

    info_file.write(out_line + "\t" + "\t".join(var_info[10:] + count_info) + "\n")


def process_format_fields(line: Sequence[str]) -> Dict[str, str]:
    """Process FORMAT fields from VCF line.

    Args:
        line: List of VCF fields

    Returns:
        Dictionary mapping FORMAT field names to their values,
        with PS field removed if present
    """
    format_values = dict(zip(line[8].split(":"), line[9].split(":")))
    if "PS" in format_values:
        del format_values["PS"]
    return format_values


def format_unassigned_variant(
    line: Sequence[str], format_values: Dict[str, str]
) -> str:
    """Format an unassigned variant line.

    Args:
        line: List of VCF fields
        format_values: Dictionary of FORMAT field values

    Returns:
        Formatted variant line string with unphased genotypes
    """
    return (
        "\t".join(
            line[0:8]
            + [":".join(format_values.keys())]
            + [":".join(format_values.values())]
        )
        .replace("1|0", "0/1")
        .replace("2|1", "1/2")
        .replace("|", "/")
    )
