from typing import Optional, Set, Tuple, Union

import pysam


def get_chroms_from_bam(filename: str) -> Set[str]:
    """Get set of chromosome names from a BAM file.

    Args:
        filename: Path to indexed BAM file

    Returns:
        Set of chromosome names from the BAM index, excluding '*' and empty entries

    Raises:
        ValueError: If BAM file is not indexed
    """
    chroms: Set[str] = set()
    stats = pysam.idxstats(filename)
    for row in stats.split("\n"):
        fields = row.split("\t")
        if fields[0] != "*" and fields[0] != "":
            chroms.add(fields[0])
    return chroms


def open_alignment(
    alignment_file: str, window: Optional[str]
) -> Tuple[
    Union[pysam.libcalignmentfile.IteratorRowRegion, str], pysam.AlignmentFile, int
]:
    """Open an alignment file and create BAM iterator for specified region.

    Args:
        alignment_file: Path to BAM/CRAM file
        window: Optional region string in format "chrom[:start[-end]]"
               Examples: "chr1", "chr1:1000", "chr1:1000-2000"

    Returns:
        Tuple containing:
            - BAM iterator for specified region (or empty string if invalid)
            - Open BAM file handle
            - Number of reads in region (0 if unspecified or invalid)

    Notes:
        Window string format supports three variants:
        - Chromosome only: "chr1"
        - Chromosome and start: "chr1:1000"
        - Chromosome, start and end: "chr1:1000-2000"
    """
    bam = pysam.AlignmentFile(alignment_file, "rb")

    if window is not None:
        window_chrom = window.split(":")[0]
        if len(window.split(":")) == 2:
            window_margin = window.split(":")[1].split("-")
            if len(window_margin) == 2:
                window_start = int(window_margin[0])
                window_end = int(window_margin[1])
                bamiter = bam.fetch(window_chrom, window_start, window_end)
                count = bam.count(window_chrom, window_start, window_end)
            else:
                window_start = int(window_margin[0])
                bamiter = bam.fetch(window_chrom, window_start)
                count = bam.count(window_chrom, window_start)
        else:
            try:
                bamiter = bam.fetch(window_chrom)
                count = bam.count(window_chrom)
            except:
                count = 0
                bamiter = ""
    else:
        bamiter = bam.fetch(until_eof=True)
        count = 0

    return bamiter, bam, count
