from patmat.io.file_utils import openfile


def get_chroms_from_vcf(vcf):
    chroms = dict()
    with openfile(vcf) as vf_file:
        for line in vf_file:
            if line.startswith("#"):
                continue
            else:
                line = line.rstrip().split("\t")
                chroms[line[0]] = int(line[1])
    return chroms
