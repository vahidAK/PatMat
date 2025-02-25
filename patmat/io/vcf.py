from typing import Dict, TextIO, Union

from patmat.io.file_utils import openfile


def get_chroms_from_vcf(vcf: str) -> Dict[str, int]:
    """Extract chromosome names and maximum positions from a VCF file.

    Parses a VCF file to create mapping of chromosome names to their
    maximum observed position.

    Args:
        vcf: Path to VCF file (can be uncompressed, .gz, or .bz2)

    Returns:
        Dictionary mapping chromosome names (str) to maximum position (int)

    Notes:
        - Skips header lines (starting with #)
        - Uses the position field (column 2) from each VCF record
        - If multiple records exist for a chromosome, keeps the largest position

    Examples:
        >>> chroms = get_chroms_from_vcf('variants.vcf.gz')
        >>> print(chroms['chr1'])  # Prints maximum position on chr1
        249250621
    """
    chroms: Dict[str, int] = dict()
    with openfile(vcf) as vf_file:
        for line in vf_file:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            chroms[line[0]] = int(line[1])
    return chroms
