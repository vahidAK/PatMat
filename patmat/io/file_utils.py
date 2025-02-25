import bz2
import gzip
from typing import TextIO, Union


def openfile(file: str) -> Union[TextIO, gzip.GzipFile, bz2.BZ2File]:
    """Open a file with appropriate compression handling.

    Opens a file for text reading, automatically handling gzip and bzip2
    compression based on file extension.

    Args:
        file: Path to file to open. Can be uncompressed, .gz, .bz or .bz2

    Returns:
        File handle appropriate for the file type:
            - TextIO for uncompressed files
            - GzipFile for .gz files
            - BZ2File for .bz/.bz2 files

    Notes:
        All files are opened in text mode ('rt')

    Examples:
        >>> f = openfile('data.txt')  # Returns TextIO
        >>> f = openfile('data.txt.gz')  # Returns GzipFile
        >>> f = openfile('data.txt.bz2')  # Returns BZ2File
    """
    if file.endswith(".gz"):
        opened_file = gzip.open(file, "rt")
    elif file.endswith("bz") or file.endswith("bz2"):
        opened_file = bz2.open(file, "rt")
    else:
        opened_file = open(file, "rt")
    return opened_file
