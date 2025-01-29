from collections import defaultdict


def read_info_iter(read_info_file):
    with open(read_info_file) as VarReadInfo:
        for line in VarReadInfo:
            line = line.rstrip().split("\t")
            for read_info in line[3].split(","):
                read_info = read_info.split(":")
                yield (line, read_info, (line[0], read_info[0]))


def build_read_dict_HP_temp(read_info_file):
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


def build_variant_dict_HP(reads_hap, read_info_file, per_var_info):
    variant_dict_HP = defaultdict(lambda: defaultdict(int))
    for line, read_info, read_key in read_info_iter(read_info_file):
        var_key = (line[0], line[1])
        hap = reads_hap.get(read_key)
        if hap:
            per_var_info[var_key][f"h{hap}all"] += 1
            if line[2] != "noninfo":
                variant_dict_HP[var_key][(line[2], hap)] += 1
    return variant_dict_HP


def add_reads_hap(hapRatio, minvariant, reads_hap, read_dict_HP_temp_reass):
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
    hap_ratio,
    min_variant,
    read_info_file,
    read_dict_HP_temp,
):

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


def build_reads_hap_temp(hapRatio, minvariant, read_dict_HP_temp):
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
