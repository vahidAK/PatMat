import pysam


def getChromsFromBAM(filename):
    chroms = set()
    stats = pysam.idxstats(filename)
    for row in stats.split("\n"):
        fields = row.split("\t")
        if fields[0] != "*" and fields[0] != "":
            chroms.add(fields[0])
    return chroms


def openalignment(alignment_file, window):
    """
    Opens an alignment file and creates bam iterator
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
