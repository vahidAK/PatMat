#! /usr/bin/env python3
# coding=utf-8

# Copyright (C) 2022  Vahid Akbari

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
PatMat: Parent-of-origin resolved chromosome-scale haplotyping
"""

__author__ = "Vahid Akbari"
__email__ = "vakbari@bcgsc.ca"
__copyright__ = "Copyright (C) 2022, " + __author__
__license__ = "GPLv3"

import os
import gzip
import bz2
import argparse
import warnings
import multiprocessing as mp
import sys
from collections import defaultdict
from itertools import repeat
import pysam
import tabix
from tqdm import tqdm

def get_variant_info(feed_list,
                  alignment_file,
                  chrom):
    read_HP_list= list()
    samfile = pysam.AlignmentFile(alignment_file, 'rb')
    for info in feed_list:
        HP,position,ref,alt= info
        try:
            sam_pileup= samfile.pileup(chrom,
                                       position,
                                       position+1,
                                       truncate=True)
        except:#The cordiniate is not found so ignore it
            warnings.warn("Variant {} {} did not find or do not have any map "
                          "reads in the alignment file. Check if correct"
                          "bam file is given or bam is indexed and not corrupted."
                          " Skipping it.".format(chrom,
                                                   position+1))
            continue
        for pileupcolumn in sam_pileup:
            pileupcolumn.set_min_base_quality(0)
            if pileupcolumn.pos == position:
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    read_id= pileupread.alignment.query_name
                    read_start= pileupread.alignment.reference_start
                    read_end= pileupread.alignment.reference_end
                    read_mq= pileupread.alignment.mapping_quality
                    read_id= pileupread.alignment.query_name
                    # read_len= pileupread.alignment.query_alignment_length
                    read_start= pileupread.alignment.reference_start
                    read_end= pileupread.alignment.reference_end
                    flag= pileupread.alignment.flag
                    phred= pileupread.alignment.query_qualities[
                                                pileupread.query_position]
                    
                    if pileupread.alignment.is_supplementary:
                        suppl_flag= str(flag)+":yes"
                    else:
                        suppl_flag= str(flag)+":no"
                    if pileupread.alignment.is_reverse:
                        strand = "-"
                    else:
                        strand = "+"
                    key_per_read = (chrom,read_start,read_end,
                                    read_id,strand,suppl_flag,read_mq)
                    val= [position,phred]
                    read_base= None
                    read_base1= None
                    read_base2= None
                    if HP == '1|0' or HP == '0|1' or HP == '0/1':
                        if len(ref) == 1 and len(alt)==1: # dealing with mismatches
                            read_base= pileupread.alignment.query_sequence[
                                    pileupread.query_position]
                        elif len(ref) > 1 and len(alt) == 1: #dealing with deletions
                            if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) -1:
                                read_base= pileupread.alignment.query_sequence[
                                        pileupread.query_position]
                            elif pileupread.indel == 0:
                                read_base = pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(ref)]
                            else:
                                continue
                        elif len(ref) == 1 and len(alt) > 1: #dealing with insertions
                            if pileupread.indel > 0 and pileupread.indel == len(alt) -1:
                                read_base= pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(alt)]
                            elif pileupread.indel == 0:
                                read_base= pileupread.alignment.query_sequence[
                                        pileupread.query_position]
                            else:
                                continue
                    elif HP == '1/2':
                        alt1,alt2 = alt
                        if len(ref) > 1 and len(alt1) > 1 and len(alt2)  > 1:
                            continue
                        if len(ref) == 1 and len(alt1) == 1: # dealing with mismatches
                            read_base1= pileupread.alignment.query_sequence[
                                    pileupread.query_position]
                        elif len(ref) > 1 and len(ref) > len(alt1): #dealing with deletions
                            if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - len(alt1):
                                read_base1= pileupread.alignment.query_sequence[
                                        pileupread.query_position:pileupread.query_position+len(ref)-abs(pileupread.indel)]
                            elif pileupread.indel == 0:
                                read_base1 = pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(ref)]
                        elif len(alt1) > 1 and len(alt1) > len(ref): #dealing with insertions
                            if pileupread.indel > 0 and pileupread.indel == len(alt1) - len(ref):
                                read_base1= pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(alt1)]
                            elif pileupread.indel == 0:
                                read_base1= pileupread.alignment.query_sequence[
                                        pileupread.query_position:pileupread.query_position+len(ref)]
                        if len(ref) == 1 and len(alt2) == 1: # dealing with mismatches
                            read_base2= pileupread.alignment.query_sequence[
                                    pileupread.query_position]
                        elif len(ref) > 1 and len(ref) > len(alt2): #dealing with deletions
                            if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) - len(alt2):
                                read_base2= pileupread.alignment.query_sequence[
                                        pileupread.query_position:pileupread.query_position+len(ref)-abs(pileupread.indel)]
                            elif pileupread.indel == 0:
                                read_base2 = pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(ref)]
                        elif len(alt2) > 1 and len(alt2) > len(ref): #dealing with insertions
                            if pileupread.indel > 0 and pileupread.indel == len(alt2) - len(ref):
                                read_base2= pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(alt2)]
                            elif pileupread.indel == 0:
                                read_base2= pileupread.alignment.query_sequence[
                                        pileupread.query_position:pileupread.query_position+len(ref)]
                                
                    if read_base is not None:
                        read_base= read_base.upper()
                        if HP == '1|0' and read_base == alt:
                            read_HP_list.append([(*key_per_read,"HP1"),
                                        ':'.join(map(str,val + [alt]))])
                        elif HP == '1|0' and read_base == ref:
                            read_HP_list.append([(*key_per_read,"HP2"),
                                        ':'.join(map(str,val + [ref]))])
                        elif HP == '0|1' and read_base == alt:
                            read_HP_list.append([(*key_per_read,"HP2"),
                                        ':'.join(map(str,val + [alt]))])
                        elif HP == '0|1' and read_base == ref:
                            read_HP_list.append([(*key_per_read,"HP1"),
                                        ':'.join(map(str,val + [ref]))])
                        elif HP == '0/1' and read_base == ref:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [ref]))])
                        elif HP == '0/1' and read_base == alt:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [alt]))])
                    if read_base1 is not None or read_base2 is not None:
                        if read_base1 is not None:
                            read_base1= read_base1.upper()
                        if read_base2 is not None:
                            read_base2= read_base2.upper()
                        alt1,alt2 = alt
                        if HP == '1/2' and read_base1 == alt1:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [alt1]))])
                        elif HP == '1/2' and read_base1 == alt2:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [alt2]))])
                        elif HP == '1/2' and read_base2 == alt1:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [alt1]))])
                        elif HP == '1/2' and read_base2 == alt2:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [alt2]))])
    return read_HP_list



def openalignment(alignment_file,
                  window):
    '''
    Opens a bam/sam file and creates bam iterator
    '''
    bam = pysam.AlignmentFile(alignment_file, 'rb')
    if window is not None:
        window_chrom = window.split(':')[0]
        if len(window.split(':')) == 2:
            window_margin= window.split(':')[1].split('-')
            if len(window_margin) == 2:
                window_start = int(window_margin[0])
                window_end = int(window_margin[1])
                bamiter = bam.fetch(window_chrom, window_start, window_end)
                count= bam.count(window_chrom, window_start, window_end)
            else:
                window_start = int(window_margin[0])
                bamiter = bam.fetch(window_chrom, window_start)
                count= bam.count(window_chrom, window_start)
        else:
            try:
                bamiter = bam.fetch(window_chrom)
                count= bam.count(window_chrom)
            except:
                count= 0
                bamiter= ""
    else:
        bamiter = bam.fetch(until_eof=True)
        count = 0
    return bamiter, bam, count


def openfile(file):
    '''
    Opens a file
    '''
    if file.endswith('.gz'):
        opened_file = gzip.open(file,'rt')
    elif file.endswith('bz') or file.endswith('bz2'):
        opened_file = bz2.open(file,'rt')
    else:
        opened_file = open(file,'rt')
    return opened_file


def PofO_dmr(known_dmr, 
             chrom_list, 
             tb_methylcall, 
             read_dictHP1,
             read_dictHP2,
             methyl_coverage,
             output,
             min_cg,
             meth_difference,
             cpg_difference):
    out_meth= open(output,'w')
    dmr_file= openfile(known_dmr)
    header= next(dmr_file).rstrip()
    out_meth.write(header+"\tAll_CpGs_InBoth_HP1-HP2_ThatMeet_CoverageCuttOff\tDifferentiallyMethylatedCpGs_HP1\t"
                   "DifferentiallyMethylatedCpGs_HP2\tMethylationFrequency_HP1\tMethylationFrequency_HP2\t"
                   "DetectionValue_HP1\tDetectionValue_HP2\n")
    chrom_hp_origin_count= defaultdict(lambda: defaultdict(int))
    for line in dmr_file:
        hp1_freq_cg = dict()
        hp2_freq_cg = dict()
        cg_sites_hp1= defaultdict(lambda: defaultdict(int))
        cg_sites_hp2= defaultdict(lambda: defaultdict(int))
        line= line.rstrip().split('\t')
        dmr_chrom= line[0]
        if dmr_chrom not in chrom_list:
            continue
        dmr_start= int(line[1]) - 1
        dmr_end= int(line[2])
        origin= line[3].lower()
        if origin not in ["maternal","paternal"]:
            warnings.warn("iDMR: {} does not have a valid origin (must be paternal or maternal. case insensitive). "
                          "This iDMR will not be used for PofO assignment and confidence calculation.".format('\t'.join(line)))
        try:
            records = tb_methylcall.query(dmr_chrom, dmr_start, dmr_end)
        except:
            warnings.warn("iDMR {}:{}-{} was ignored because it does not have any reads in "
                          "methylation call file.".format(dmr_chrom, dmr_start, dmr_end))
            records = "NA"
        if records != "NA":
            for record in records:
                key= (record[0],record[4],record[3])
                if key in read_dictHP1[dmr_chrom] and key not in read_dictHP2[dmr_chrom]:
                    if record[7] != "NA":
                        for mod_sites in record[7].split(','):
                            mod_sites= int(mod_sites)
                            if mod_sites >= dmr_start and mod_sites <= dmr_end:
                                cg_sites_hp1["mod"][mod_sites] += 1
                                cg_sites_hp1["all"][mod_sites] += 1
                    if record[8] != "NA":
                        for unmod_sites in record[8].split(','):
                            unmod_sites = int(unmod_sites)
                            if unmod_sites >= dmr_start and unmod_sites <= dmr_end:
                                cg_sites_hp1["all"][unmod_sites] += 1
                elif key in read_dictHP2[dmr_chrom] and key not in read_dictHP1[dmr_chrom]:
                    if record[7] != "NA":
                        for mod_sites in record[7].split(','):
                            mod_sites= int(mod_sites)
                            if mod_sites >= dmr_start and mod_sites <= dmr_end:
                                cg_sites_hp2["mod"][mod_sites] += 1
                                cg_sites_hp2["all"][mod_sites] += 1
                    if record[8] != "NA":
                        for unmod_sites in record[8].split(','):
                            unmod_sites = int(unmod_sites)
                            if unmod_sites >= dmr_start and unmod_sites <= dmr_end:
                                cg_sites_hp2["all"][unmod_sites] += 1
            for cg_site,num_all in cg_sites_hp1['all'].items():
                try:
                    num_mod = cg_sites_hp1['mod'][cg_site]
                except:
                    num_mod = 0
                if num_all >= methyl_coverage:
                    hp1_freq_cg[cg_site]= num_mod/num_all
                
            for cg_site,num_all in cg_sites_hp2['all'].items():
                try:
                    num_mod = cg_sites_hp2['mod'][cg_site]
                except:
                    num_mod = 0
                if num_all >= methyl_coverage:
                    hp2_freq_cg[cg_site]= num_mod/num_all
        common_cg= [i for i in list(hp1_freq_cg.keys()) if i in list(hp2_freq_cg.keys())]
        diff_cg_hp1= len([i for i in common_cg if hp1_freq_cg[i] - hp2_freq_cg[i] >= meth_difference])
        diff_cg_hp2= len([i for i in common_cg if hp2_freq_cg[i] - hp1_freq_cg[i] >= meth_difference])
        hp1_freq= 0
        hp2_freq= 0
        for i in common_cg:
            hp1_freq += hp1_freq_cg[i] 
            hp2_freq += hp2_freq_cg[i] 
        if len(common_cg) < 1:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+"NA"+'\t'+"NA"+
                           '\t'+"NA" + '\t' + "NA"+'\t' + "Ignored:No_CpG" + '\t' +
                           "Ignored:No_CpG" +'\n')
            continue
        hp1_freq= hp1_freq/len(common_cg)
        hp2_freq= hp2_freq/len(common_cg)
        if len(common_cg) < min_cg and abs(diff_cg_hp1 - diff_cg_hp2) / len(common_cg) < cpg_difference:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg_And_cpg_difference" + '\t' +
                           "Ignored:DidNotMeet_min_cg_And_cpg_difference" +'\n')
            continue
        elif len(common_cg) < min_cg and diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" + '\t' +
                           "Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" +'\n')
            continue
        elif len(common_cg) < min_cg:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg" + '\t' +
                           "Ignored:DidNotMeet_min_cg" +'\n')
            continue
        elif abs(diff_cg_hp1 - diff_cg_hp2) / len(common_cg) < cpg_difference:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_cpg_difference" + '\t' +
                           "Ignored:DidNotMeet_cpg_difference" +'\n')
            continue
        elif diff_cg_hp1 == 0 and diff_cg_hp2 == 0: 
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+ 
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" + '\t' +
                           "Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" +'\n') 
            continue

        num_cg_to_add_hp1= (diff_cg_hp1 * hp1_freq)/len(common_cg)
        num_cg_to_add_hp2= (diff_cg_hp2 * hp2_freq)/len(common_cg)
        out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                       '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + str(num_cg_to_add_hp1) + '\t' + 
                       str(num_cg_to_add_hp2)+'\n')  
        if origin == 'maternal':
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += num_cg_to_add_hp1
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += num_cg_to_add_hp2
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += num_cg_to_add_hp1
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += num_cg_to_add_hp2
        if origin == 'paternal':
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += num_cg_to_add_hp2 
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += num_cg_to_add_hp1
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += num_cg_to_add_hp2 
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += num_cg_to_add_hp1
    dmr_file.close()
    out_meth.close()
    return chrom_hp_origin_count

def vcf2dict(vcf_strand):
    """
    This function converts the input vcf file to haplotype1 and
    haplotype2 dictionaries to be used for read phasing.
    """
    vcf_dict= defaultdict(dict)
    vcf_file = openfile(vcf_strand)
    for line in vcf_file:
        if line.startswith("#"):
            continue
        line_list = line.rstrip().split('\t')
        if (line_list[9].startswith('.|0') or 
            line_list[9].startswith('0|.') or
            line_list[9].startswith('1|.') or
            line_list[9].startswith('.|1')):
            line_list[9] = line_list[9].replace(".|0","1|0").replace(".|1","0|1").replace("1|.","1|0").replace("0|.","0|1")
            warnings.warn("{}:{} variant in strand-seq vcf has .|0 or 0|. or .|1 or 1|. phased genotype."
                          "Note that it will be interpreted as 1|0 or 0|1 or 0|1 or 1|0".format(line_list[0], line_list[1]))
        chrom = line_list[0]
        if (line_list[9].startswith('1|0') or 
            line_list[9].startswith('0|1')):
            vcf_dict[chrom][line_list[1]] = line_list[9].split(':')[0]
    vcf_file.close()
    return vcf_dict


def strand_vcf2dict_phased(vcf_strand, vcf):
    final_dict= defaultdict(set)
    vcf_dict= vcf2dict(vcf_strand)
    with openfile(vcf) as vf:
        for line in vf:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            if ((line[9].startswith('0/1') or 
                 line[9].startswith('1/0') or 
                 line[9].startswith('0|1') or
                 line[9].startswith('1|0'))
                and
                (line[0] in vcf_dict and 
                 line[1] in vcf_dict[line[0]])):
                final_dict[line[0]].add((vcf_dict[line[0]][line[1]],
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif (line[9].startswith('0/1') or 
                  line[9].startswith('1/0') or 
                  line[9].startswith('0|1') or
                  line[9].startswith('1|0')):
                final_dict[line[0]].add(("0/1",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif (line[9].startswith('1/2') or 
                  line[9].startswith('1|2') or
                  line[9].startswith('2/1') or 
                  line[9].startswith('2|1')):
                final_dict[line[0]].add(("1/2",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        (line[4].split(',')[0].upper(),
                                        line[4].split(',')[1].upper())))
    return final_dict
                

def vcf2dict_phased(block_file, 
                    vcf_whats, 
                    vcf_strand,
                    hapRatio,
                    minvariant,
                    vcf):
    """
    This function converts the input vcf file to haplotype1 and
    haplotype2 dictionaries to be used for read phasing.
    """
    final_dict= defaultdict(set)
    vcf_dict= vcf2dict(vcf_strand)
    tb_whatsvcf= tabix.open(os.path.abspath(vcf_whats))
    for blocks in block_file:
        b_chrom, b_start, b_end= blocks
        agreement_count= 0
        disagreement_count= 0
        try:
            records_whats = tb_whatsvcf.query(b_chrom, b_start-1, b_end+1)
            records_whats = list(records_whats)
        except:
            warnings.warn("{}:{}-{} block cannot be extracted from WhatsHap vcf file. "
                          "Make sure file is indexed. Skipping it.".format(b_chrom, b_start, b_end))
            continue
        if b_chrom in vcf_dict:
            disagreement_sites = set()
            agreement_sites = set()
            for vcf_line in records_whats:
                if vcf_line[1] in vcf_dict[vcf_line[0]]:
                    if vcf_line[9].startswith("1|0") and vcf_dict[vcf_line[0]][vcf_line[1]] == "1|0":
                        agreement_sites.add(vcf_line[1])
                        agreement_count += 1
                    elif vcf_line[9].startswith("1|0") and vcf_dict[vcf_line[0]][vcf_line[1]] == "0|1":
                        disagreement_sites.add(vcf_line[1])
                        disagreement_count += 1
                    elif vcf_line[9].startswith("0|1") and vcf_dict[vcf_line[0]][vcf_line[1]] == "0|1":
                        agreement_sites.add(vcf_line[1])
                        agreement_count += 1
                    elif vcf_line[9].startswith("0|1") and vcf_dict[vcf_line[0]][vcf_line[1]] == "1|0":
                        disagreement_sites.add(vcf_line[1])
                        disagreement_count += 1                    
            for vcf_line in records_whats:
                if vcf_line[1] in vcf_dict[vcf_line[0]]:
                    continue
                if (vcf_line[9].startswith('1|0') or 
                    vcf_line[9].startswith('0|1')):
                    if (agreement_count > disagreement_count and
                        agreement_count >= minvariant and 
                        agreement_count/(agreement_count+disagreement_count) >= hapRatio and
                        vcf_line[1] not in disagreement_sites):
                        vcf_dict[vcf_line[0]][vcf_line[1]] = vcf_line[9].split(':')[0]
                    elif (disagreement_count > agreement_count and
                          disagreement_count >= minvariant and 
                          disagreement_count/(agreement_count+disagreement_count) >= hapRatio and
                          vcf_line[1] not in agreement_sites):
                        vcf_dict[vcf_line[0]][vcf_line[1]] = vcf_line[9].split(':')[0][::-1]

    with openfile(vcf) as vf:
        for line in vf:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            if ((line[9].startswith('0/1') or 
                 line[9].startswith('1/0') or 
                 line[9].startswith('0|1') or
                 line[9].startswith('1|0'))
                and
                (line[0] in vcf_dict and 
                 line[1] in vcf_dict[line[0]])):
                final_dict[line[0]].add((vcf_dict[line[0]][line[1]],
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif (line[9].startswith('0/1') or 
                  line[9].startswith('1/0') or 
                  line[9].startswith('0|1') or
                  line[9].startswith('1|0')):
                final_dict[line[0]].add(("0/1",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif (line[9].startswith('1/2') or 
                  line[9].startswith('1|2') or
                  line[9].startswith('2/1') or 
                  line[9].startswith('2|1')):
                final_dict[line[0]].add(("1/2",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        (line[4].split(',')[0].upper(),
                                        line[4].split(',')[1].upper())))

    vcf_dict.clear()
    return final_dict


def per_read_variant(vcf_dict,
                     bam_file,
                     chunk,
                     threads,
                     perReadinfo):
    chrom_list = sorted(list(vcf_dict.keys()))
    for chrom in chrom_list:
        per_read_hp = defaultdict(lambda: defaultdict(list))
        bamiter, bam, count = openalignment(bam_file, chrom)
        if count > 0:
            feed_list= list(vcf_dict[chrom])
            feed_list = [feed_list[x:x+chunk]
                                 for x in range(0, len(feed_list),
                                                chunk)]
            feed_list = [feed_list[x:x+threads]
                                 for x in range(0, len(feed_list),
                                                threads)]
            description= "Tagging variants to reads from {}: ".format(chrom)
            with tqdm(total=len(feed_list),
                desc=description,
                bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                                  ) as pbar:
                for vcf_info_list in feed_list:
                    p= mp.Pool(len(vcf_info_list))
                    results= p.starmap(get_variant_info,
                                       list(zip(vcf_info_list,
                                                repeat(bam_file),
                                                repeat(chrom))))
                    p.close()
                    p.join()
                    for result in results:
                        if result is not None:
                            for read_info in result:
                                key,val= read_info
                                per_read_hp[key[0:-1]][key[-1]].append(val)
                    pbar.update(1)
            for key in list(per_read_hp.keys()):
                try:
                    hp1_variants= per_read_hp[key]['HP1']
                    hp1_count= str(len(hp1_variants))
                except:
                    hp1_count= '0'
                if hp1_count == '0':
                    hp1_variants= ['NA']
                    
                try:
                    hp2_variants= per_read_hp[key]['HP2']
                    hp2_count= str(len(hp2_variants))
                except:
                    hp2_count= '0'   
                if hp2_count == '0':
                    hp2_variants= ['NA']
                    
                try:
                    unphased_variants= per_read_hp[key]['unphased']
                    unphased_count= str(len(unphased_variants))
                except:
                    unphased_count= '0'
                if unphased_count == '0':
                    unphased_variants= ['NA']
                    
                out_to_write= list(map(str,key)) + [','.join(sorted(hp1_variants)), 
                                                    ','.join(sorted(hp2_variants)),
                                                    ','.join(sorted(unphased_variants))]
                
                perReadinfo.write('\t'.join(out_to_write)+'\n')
        else:
            warnings.warn("{} does not have any mapped reads in alignment "
                          "file Or alignment is truncated or corrupt indexed. "
                          "Skipping it.".format(chrom))
        per_read_hp.clear()


def get_block(vcf_whats):
    blocks= defaultdict(set)
    final= list()
    with openfile(vcf_whats) as wv:
        for line in wv:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            if "|" in line[9].split(":")[0]:
               blocks[(line[0],line[9].split(":")[-1])].add(int(line[1]))
    for key, val in blocks.items():
        val= sorted(val)
        final.append((key[0], val[0], val[-1]))
    blocks.clear()
    return final


def get_indels(vcf):
    indels= set()
    with openfile(vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue
            line= line.rstrip().split('\t')
            if (line[9].startswith("0|1") or
                line[9].startswith("1|0") or
                line[9].startswith("0/1") or
                line[9].startswith("1/0")):
                if len(line[3]) > 1 or len(line[4]) > 1:
                    indels.add((line[0],str(int(line[1])-1)))
            # Non-reference het indels (e.g. 1/2) are ignored as they will be in unphased column of per-read info
    return indels



def main(args):
    """
    This is the phase module which phase the nanopore reads and
    methylation data to corresponding haplotype using vcf file and
    processed methylation call file.
    """
    hapRatio = args.hapratio
    minvariant= args.min_variant
    min_cg= args.min_cg
    bam_file = os.path.abspath(args.bam)
    threads = args.threads
    methyl_coverage= args.methyl_coverage
    meth_difference= args.meth_difference
    cpg_difference= args.cpg_difference
    MinBaseQuality = args.min_base_quality
    chunk = args.chunk_size
    min_read_reassignment= args.min_read_number
    out = os.path.abspath(args.output)
    MappingQuality = args.mapping_quality
    vcf= os.path.abspath(args.vcf)
    known_dmr= os.path.abspath(args.known_dmr)
    MethylCallfile = os.path.abspath(args.methylcallfile)
    tb_methylcall = tabix.open(MethylCallfile)
                
    if not os.path.isfile(MethylCallfile+".tbi"):
        raise Exception("It seems processed methylation call file "
                        "is not index by tabix.")
    sites_to_ignore= set()
    if args.black_list is not None:
        if not os.path.isfile(vcf+".tbi"):
            raise Exception("Black list is given but it seems that"
                            " the vcf file is not indexed (file ends with"
                            " .tbi was not found).")
        tb_vcf= tabix.open(vcf)
        black_list= os.path.abspath(args.black_list)
        with openfile(black_list) as bl:
            for line in bl:
                line= line.rstrip().split('\t')
                try:
                    records= tb_vcf.query(line[0], int(line[1]), int(line[2])+1)
                except:
                    warnings.warn("{}:{}-{} region from black list does not exist in the "
                          "vcf file. Skipping it.".format(line[0],line[1] , line[2]))
                    records= "NA"
                if records != "NA":
                    for record in records:
                        sites_to_ignore.add((record[0],str(int(record[1])-1)))
    if args.per_read is not None:
        per_read_file= os.path.abspath(args.per_read)
    else:
        if args.strand_vcf is not None and args.whatshap_vcf is None:
            warnings.warn("Using strand-seq phased vcf only with {}.".format(known_dmr.split('/')[-1]))
            vcf_strand = os.path.abspath(args.strand_vcf)
            final_dict= strand_vcf2dict_phased(vcf_strand, vcf)
        # elif args.strand_vcf is None and args.whatshap_vcf is not None:
        #     warnings.warn("No strand-seq phased vcf is given. Using WhatsHap phased vcf only.")
        #     vcf_whats= os.path.abspath(args.whatshap_vcf)
        #     if args.whatshap_block is None:
        #         block_file= get_block(vcf_whats)
        #     else:
        #         block_file= list()
        #         with openfile(args.whatshap_block) as wb:
        #             next(wb)
        #             for line in wb:
        #                 line=line.rstrip().split('\t')
        #                 block_file.append((line[0],int(line[1]),int(line[2])))
        elif args.strand_vcf is not None and args.whatshap_vcf is not None:
            warnings.warn("Using both strand-seq phased and WhatsHap phased vcf with {}.".format(known_dmr.split('/')[-1]))
            vcf_whats= os.path.abspath(args.whatshap_vcf)
            vcf_strand = os.path.abspath(args.strand_vcf)
            if not os.path.isfile(vcf_whats+".tbi"):
                raise Exception("It seems that whatshap vcf "
                                "is not index by tabix.")
            if args.whatshap_block is None:
                block_file= get_block(vcf_whats)
            else:
                block_file= list()
                with openfile(args.whatshap_block) as wb:
                    next(wb)
                    for line in wb:
                        line=line.rstrip().split('\t')
                        block_file.append((line[0],int(line[1]),int(line[2])))
            final_dict= vcf2dict_phased(block_file, 
                                        vcf_whats, 
                                        vcf_strand,
                                        hapRatio,
                                        minvariant,
                                        vcf)
        else:
            raise Exception("No strand-seq vcf is given.")
    if args.per_read is None:
        out_per_read = out + '_HP1'
        perReadinfo= open(out_per_read+"_HP2_PerReadInfo.tsv", 'w')
        perReadinfo.write("#Chromosome\tReadRefStart\tReadRefEnd\tReadID\t"
                          "Strand\tReadFlag:Is_Supplementary\tReadMapQuality\t"
                          "Position:BaseQuality:HP1-variants\t"
                          "Position:BaseQuality:HP2-variants\t"
                          "Position:BaseQuality:UnPhasedAndOtherHetvariants\n")
        per_read_variant(final_dict,
                         bam_file,
                         chunk,
                         threads,
                         perReadinfo)
        final_dict.clear()
        perReadinfo.close()
        per_read= openfile(out_per_read+"_HP2_PerReadInfo.tsv")
    else:
        per_read= openfile(per_read_file)
    chrom_list= set()
    read_dictHP1 = defaultdict(set)
    read_dictHP2 = defaultdict(set)
    ignore_indels= set()
    if not args.include_indels:
        ignore_indels= get_indels(vcf)
    for line in per_read:
        if line.startswith("#"):
            continue
        line= line.rstrip().split('\t')
        if int(line[6]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        chrom_list.add(line[0])
        key= (line[0],line[3],line[4])
        hp1s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[7].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality and 
                      (line[0],i.split(":")[0]) not in sites_to_ignore and
                      (line[0],i.split(":")[0]) not in ignore_indels]
        hp2s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[8].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality and 
                      (line[0],i.split(":")[0]) not in sites_to_ignore and
                      (line[0],i.split(":")[0]) not in ignore_indels]
        hp1_count= len(hp1s)
        hp2_count= len(hp2s)
        if (hp1_count > hp2_count and 
            hp1_count/(hp1_count+hp2_count) >= hapRatio and 
            hp1_count >= minvariant):
            read_dictHP1[line[0]].add(key)
            
        elif (hp2_count > hp1_count and 
              hp2_count/(hp1_count+hp2_count) >= hapRatio and 
              hp2_count >= minvariant):
            read_dictHP2[line[0]].add(key)
    per_read.close()
    variant_dict_HP1= defaultdict(lambda: defaultdict(int))
    variant_dict_HP2= defaultdict(lambda: defaultdict(int))
    if args.per_read is not None:
        per_read= openfile(os.path.abspath(args.per_read))
    else:
        per_read= openfile(out_per_read+"_HP2_PerReadInfo.tsv")
    for line in per_read:
        if line.startswith("#"):
            continue
        line= line.rstrip().split('\t')
        if int(line[6]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        key= (line[0],line[3],line[4])
        variants= [(line[0],i.split(":")[0],i.split(':')[2]) for i in 
                    set(line[7].split(',') + line[8].split(',') +line[9].split(','))
                    if i != 'NA' and int(i.split(':')[1]) >= MinBaseQuality]
        if key in read_dictHP1[line[0]]:
            for i in variants:
                variant_dict_HP1[(i[0],i[1])][i[2]] += 1
                variant_dict_HP2[(i[0],i[1])][i[2]] += 0
        elif key in read_dictHP2[line[0]]:
            for i in variants:
                variant_dict_HP1[(i[0],i[1])][i[2]] += 0
                variant_dict_HP2[(i[0],i[1])][i[2]] += 1
    per_read.close()

    if not read_dictHP1 and not read_dictHP2:
        raise Exception("No reads could be phased. Probably phased vcf files"
                        " do not have any phased variants or have very few.")
    out_non_pofo = out + '_NonPofO_HP1-HP2_Reassignment.vcf'
    out_non_pofo_reads = out + '_NonPofO_HP1-HP2_reads.tsv'
    re_assignment= open(out_non_pofo,'w')
    reads_NonPofO= open(out_non_pofo_reads,'w')
    hp1s= set()
    hp2s= set()
    with openfile(vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                re_assignment.write(line)
                continue
            line= line.rstrip().split('\t')
            if 'PS' in line[8].split(":"):
                ps_index= line[8].split(":").index('PS')
                new_ps= line[8].split(":")
                new_hp= line[9].split(":")
                new_ps.pop(ps_index)
                new_hp.pop(ps_index)
            else:
                new_ps= line[8].split(":")
                new_hp= line[9].split(":")
            if line[0] in chrom_list:
                if (line[9].startswith("0/1") or
                    line[9].startswith("1/0") or
                    line[9].startswith("0|1") or
                    line[9].startswith("1|0")):
                    hp1_count_alt= variant_dict_HP1[(line[0],str(int(line[1])-1))][line[4].upper()]
                    hp2_count_alt= variant_dict_HP2[(line[0],str(int(line[1])-1))][line[4].upper()]
                    hp1_count_ref= variant_dict_HP1[(line[0],str(int(line[1])-1))][line[3].upper()]
                    hp2_count_ref= variant_dict_HP2[(line[0],str(int(line[1])-1))][line[3].upper()]
                    if ((hp1_count_alt > hp2_count_alt and
                         hp1_count_alt > hp1_count_ref and
                         hp1_count_alt >= min_read_reassignment) or
                        (hp2_count_ref > hp1_count_ref and
                         hp2_count_ref > hp2_count_alt and
                         hp2_count_ref >= min_read_reassignment)):
                        if not args.include_indels:
                            if len(line[3]) == 1 and len(line[4]) == 1:
                                hp1s.add((line[0],line[1],line[4].upper()))
                                hp2s.add((line[0],line[1],line[3].upper()))
                        else:
                            hp1s.add((line[0],line[1],line[4].upper()))
                            hp2s.add((line[0],line[1],line[3].upper()))    
                        re_assignment.write('\t'.join(line[0:8]+
                                                      [':'.join(new_ps)+":PS"]+
                                                      ["1|0:"+':'.join(new_hp[1:])+":HP2|HP1"])+'\n')
                    elif ((hp2_count_alt > hp1_count_alt and
                           hp2_count_alt > hp2_count_ref and
                           hp2_count_alt >= min_read_reassignment) or
                          (hp1_count_ref > hp2_count_ref and
                           hp1_count_ref > hp1_count_alt and
                           hp1_count_ref >= min_read_reassignment)):
                        if not args.include_indels:
                            if len(line[3]) == 1 and len(line[4]) == 1:
                                hp1s.add((line[0],line[1],line[3].upper()))
                                hp2s.add((line[0],line[1],line[4].upper()))
                        else:
                            hp1s.add((line[0],line[1],line[3].upper()))
                            hp2s.add((line[0],line[1],line[4].upper()))
                        re_assignment.write('\t'.join(line[0:8]+
                                                      [':'.join(new_ps)+":PS"]+
                                                      ["0|1:"+':'.join(new_hp[1:])+":HP1|HP2"])+'\n')
                    else:
                        re_assignment.write('\t'.join(line[0:8]+
                                                      [':'.join(new_ps)]+
                                                      [':'.join(new_hp).replace("|", "/")])+'\n')
                elif line[9].startswith("1/1") or line[9].startswith("1|1"):
                    re_assignment.write('\t'.join(line)+'\n')
                
                elif (line[9].startswith('1/2') or 
                      line[9].startswith('1|2') or
                      line[9].startswith('2/1') or 
                      line[9].startswith('2|1')):
                    hp1_count_alt= variant_dict_HP1[(line[0],str(int(line[1])-1))][line[4].split(',')[1].upper()]
                    hp2_count_alt= variant_dict_HP2[(line[0],str(int(line[1])-1))][line[4].split(',')[1].upper()]
                    hp1_count_ref= variant_dict_HP1[(line[0],str(int(line[1])-1))][line[4].split(',')[0].upper()]
                    hp2_count_ref= variant_dict_HP2[(line[0],str(int(line[1])-1))][line[4].split(',')[0].upper()]
                    if ((hp1_count_alt > hp2_count_alt and
                         hp1_count_alt > hp1_count_ref and
                         hp1_count_alt >= min_read_reassignment) or
                        (hp2_count_ref > hp1_count_ref and
                         hp2_count_ref > hp2_count_alt and
                         hp2_count_ref >= min_read_reassignment)):
                        if not args.include_indels:
                            if len(line[3]) == 1 and len(line[4]) == 3:
                                hp1s.add((line[0],line[1],line[4].split(',')[1].upper()))
                                hp2s.add((line[0],line[1],line[4].split(',')[0].upper()))
                        else:
                            hp1s.add((line[0],line[1],line[4].split(',')[1].upper()))
                            hp2s.add((line[0],line[1],line[4].split(',')[0].upper()))
                        re_assignment.write('\t'.join(line[0:8]+
                                            [':'.join(new_ps)+":PS"]+
                                            ["1|2:"+':'.join(new_hp[1:])+":Ref_HP2|HP1"])+'\n')
                    elif ((hp2_count_alt > hp1_count_alt and
                           hp2_count_alt > hp2_count_ref and
                           hp2_count_alt >= min_read_reassignment) or
                          (hp1_count_ref > hp2_count_ref and
                           hp1_count_ref > hp1_count_alt and
                           hp1_count_ref >= min_read_reassignment)):
                        if not args.include_indels:
                            if len(line[3]) == 1 and len(line[4]) == 3:
                                hp1s.add((line[0],line[1],line[4].split(',')[0].upper()))
                                hp2s.add((line[0],line[1],line[4].split(',')[1].upper()))
                        else:
                            hp1s.add((line[0],line[1],line[4].split(',')[0].upper()))
                            hp2s.add((line[0],line[1],line[4].split(',')[1].upper()))
                        re_assignment.write('\t'.join(line[0:8]+
                                            [':'.join(new_ps)+":PS"]+
                                            ["1|2:"+':'.join(new_hp[1:])+":Ref_HP1|HP2"])+'\n')
                    else:
                        re_assignment.write('\t'.join(line[0:8]+
                                            [':'.join(new_ps)]+
                                            [':'.join(new_hp).replace("|", "/")])+'\n')
        re_assignment.close()
    
    read_dictHP1 = defaultdict(set)
    read_dictHP2 = defaultdict(set)    
    if args.per_read is not None:
        per_read= openfile(os.path.abspath(args.per_read))
    else:
        per_read= openfile(out_per_read+"_HP2_PerReadInfo.tsv")
    for line in per_read:
        if line.startswith("#"):
            continue
        line= line.rstrip().split('\t')
        if int(line[6]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        key= (line[0],line[3],line[4])
        variants= [(line[0],str(int(i.split(":")[0])+1),i.split(':')[2]) for i in 
                        set(line[7].split(',') + line[8].split(',') +line[9].split(',')) 
                        if i != 'NA' and int(i.split(':')[1]) >= MinBaseQuality]
        hp1_count= len([i for i in variants if i in hp1s])
        hp2_count= len([i for i in variants if i in hp2s])
        if (hp1_count > hp2_count and 
            hp1_count/(hp1_count+hp2_count) >= hapRatio and 
            hp1_count >= minvariant):
            read_dictHP1[line[0]].add(key)
        elif (hp2_count > hp1_count and 
              hp2_count/(hp1_count+hp2_count) >= hapRatio and 
              hp2_count >= minvariant):
            read_dictHP2[line[0]].add(key)
    per_read.close()
    for chrom, reads in read_dictHP1.items():
        for read in reads:
            reads_NonPofO.write('\t'.join(read)+'\t'+"HP1"+'\n')
    for chrom, reads in read_dictHP2.items():
        for read in reads:
            reads_NonPofO.write('\t'.join(read)+'\t'+"HP2"+'\n')
    reads_NonPofO.close()

    if args.known_dmr is not None:
        out_pofo = out + '_PofO_Assignment.vcf'
        out_pofo_reads = out + '_PofO_Assignment_reads.tsv'
        assignment_file= open(out_pofo, 'w')
        reads_PofO= open(out_pofo_reads,'w')
        chrom_list= sorted(chrom_list)
        chrom_hp_origin_count= PofO_dmr(known_dmr, 
                                        chrom_list, 
                                        tb_methylcall, 
                                        read_dictHP1,
                                        read_dictHP2,
                                        methyl_coverage,
                                        out + '_CpG-Methylation-Status-at-DMRs.tsv',
                                        min_cg,
                                        meth_difference,
                                        cpg_difference)
        chrom_hp_origin = defaultdict(dict)
        if chrom_hp_origin_count:
            for key, val in chrom_hp_origin_count.items():
                chrom, origin= key
                hp1_mat_count= 0
                hp2_mat_count= 0
                hp1_pat_count= 0
                hp2_pat_count= 0
                if origin.lower() == 'maternal':
                    hp1_mat_count= val['HP1']
                    hp2_mat_count= val['HP2']
                    if hp1_mat_count > hp2_mat_count:
                        chrom_hp_origin[chrom]['HP1'] = ('maternal',
                                                          hp1_mat_count/(hp1_mat_count + hp2_mat_count))
                        chrom_hp_origin[chrom]['HP2'] = ('paternal',
                                                          hp1_mat_count/(hp1_mat_count + hp2_mat_count))
                    elif hp2_mat_count > hp1_mat_count:
                        chrom_hp_origin[chrom]['HP2'] = ('maternal',
                                                          hp2_mat_count/(hp1_mat_count + hp2_mat_count))
                        chrom_hp_origin[chrom]['HP1'] = ('paternal',
                                                          hp2_mat_count/(hp1_mat_count + hp2_mat_count))
                elif origin.lower() == 'paternal':
                    hp1_pat_count= val['HP1']
                    hp2_pat_count= val['HP2']
                    if hp1_pat_count > hp2_pat_count:
                        chrom_hp_origin[chrom]['HP1'] = ('paternal',
                                                          hp1_pat_count/(hp1_pat_count + hp2_pat_count))
                        chrom_hp_origin[chrom]['HP2'] = ('maternal',
                                                          hp1_pat_count/(hp1_pat_count + hp2_pat_count))
                    elif hp2_pat_count > hp1_pat_count:
                        chrom_hp_origin[chrom]['HP2'] = ('paternal',
                                                        hp2_pat_count/(hp1_pat_count + hp2_pat_count))
                        chrom_hp_origin[chrom]['HP1'] = ('maternal',
                                                        hp2_pat_count/(hp1_pat_count + hp2_pat_count))
        if chrom_hp_origin:
            for chrom, reads in read_dictHP1.items():
                if chrom not in chrom_hp_origin:
                    continue
                if chrom_hp_origin[chrom]['HP1'][0] == 'maternal':
                    for read in reads:
                        reads_PofO.write('\t'.join(read)+'\t'+"maternal"+'\n')
                if chrom_hp_origin[chrom]['HP1'][0] == 'paternal':
                    for read in reads:
                        reads_PofO.write('\t'.join(read)+'\t'+"paternal"+'\n')
            for chrom, reads in read_dictHP2.items():
                if chrom not in chrom_hp_origin:
                    continue
                if chrom_hp_origin[chrom]['HP2'][0] == 'maternal':
                    for read in reads:
                        reads_PofO.write('\t'.join(read)+'\t'+"maternal"+'\n')
                if chrom_hp_origin[chrom]['HP2'][0] == 'paternal':
                    for read in reads:
                        reads_PofO.write('\t'.join(read)+'\t'+"paternal"+'\n')
            reads_PofO.close()
            with openfile(out_non_pofo) as vcf_file:
                for line in vcf_file:
                    if line.startswith('#'):
                        assignment_file.write(line)
                        continue
                    line= line.rstrip().split('\t')
                    if line[0] in chrom_hp_origin:
                        if line[9].startswith('0|1'):
                            if chrom_hp_origin[line[0]]['HP1'][0] == 'maternal':
                                assignment_file.write('\t'.join(line[0:9])+'\t'+
                                                          line[9].replace("HP1","Mat").replace("HP2","Pat")
                                                          +'\n')
                            if chrom_hp_origin[line[0]]['HP1'][0] == 'paternal':
                                assignment_file.write('\t'.join(line[0:9])+'\t'+
                                                          line[9].replace("0|1","1|0").replace("HP1","Pat").replace("HP2","Mat")
                                                          +'\n')
                        elif line[9].startswith('1|0'):
                            if chrom_hp_origin[line[0]]['HP2'][0] == 'maternal':
                                assignment_file.write('\t'.join(line[0:9])+'\t'+
                                                          line[9].replace("1|0","0|1").replace("HP1","Pat").replace("HP2","Mat")
                                                          +'\n')
                            if chrom_hp_origin[line[0]]['HP2'][0] == 'paternal':
                                assignment_file.write('\t'.join(line[0:9])+'\t'+
                                                          line[9].replace("HP1","Mat").replace("HP2","Pat")
                                                          +'\n')
                        elif line[9].startswith('1|2'):
                            if chrom_hp_origin[line[0]]['HP2'][0] == 'maternal':
                                assignment_file.write('\t'.join(line[0:9])+"\t"+
                                                          line[9].replace("HP2", "Mat").replace("HP1", "Pat")
                                                          +'\n')
                            if chrom_hp_origin[line[0]]['HP2'][0] == 'paternal':
                                assignment_file.write('\t'.join(line[0:9])+'\t'+
                                                          line[9].replace("HP2", "Pat").replace("HP1", "Mat")
                                                          +'\n')
                        else:
                            assignment_file.write('\t'.join(line[0:8]+
                                                    [line[8].replace(":PS", "")]+
                                                    [line[9].replace("|", "/").
                                                     replace(":Ref_HP1/HP2","").
                                                     replace(":Ref_HP2/HP1","").
                                                     replace(":HP1/HP2","").
                                                     replace(":HP2/HP1","")])+'\n')
                    else:
                        assignment_file.write('\t'.join(line[0:8]+
                                                    [line[8].replace(":PS", "")]+
                                                    [line[9].replace("|", "/").
                                                      replace(":Ref_HP1/HP2","").
                                                      replace(":Ref_HP2/HP1","").
                                                      replace(":HP1/HP2","").
                                                      replace(":HP2/HP1","")])+'\n')
        assignment_file.close()
        print("Chromosome\tOrigin_HP1\tOrigin_HP2\tConfidence")
        for chrom,val in chrom_hp_origin.items():
            for hp,score in val.items():
                if hp=="HP1":
                    origin_hp1= score[0]
                    if origin_hp1 == "maternal":
                        origin_hp1= "Maternal"
                        origin_hp2= "Paternal"
                    elif origin_hp1 == "paternal":
                        origin_hp1= "Paternal"
                        origin_hp2= "Maternal"
                    print('\t'.join([chrom,origin_hp1,origin_hp2,str(score[1])]))
    else:
        sys.stderr.write("known DMR file and methylation call "
                                      "file are not given to determin PofO.\n")
 
    

"""
Specific argument parser.
"""
parser = argparse.ArgumentParser(prog='PatMat.py',
                                  description="Phasing reads and Methylation "
                                  "using strand-seq and nanopore to determine "
                                  "PofO of each homologous chromosome "
                                  "in a single sample.")
sp_input = parser.add_argument_group("required arguments")
sp_input.add_argument("--bam", "-b",
                      action="store",
                      type=str,
                      required=True,
                      help="The path to the cordinate sorted bam file.")
sp_input.add_argument("--output", "-o",
                      action="store",
                      type=str,
                      required=True,
                      help=("The path to directory and prefix to save "
                            "files. e.g path/to/directory/prefix"))
sp_input.add_argument("--vcf", "-v",
                  action="store",
                  type=str,
                  required=True,
                  default= None,
                  help="The path to the vcf file.")
sp_input.add_argument("--strand_vcf", "-sv",
                      action="store",
                      type=str,
                      required=True,
                      help="The path to the chromosome-scale phased vcf file."
                           "This is the input vcf file that has been phased using strand-seq data.")
sp_input.add_argument("--methylcallfile", "-mc",
                      action="store",
                      type=str,
                      required=True,
                      help=("The path to the bgziped and indexed methylation "
                            "call file processed using NanoMethPhase"
                            " methyl_call_processor Module."))
sp_input = parser.add_argument_group("Optional arguments.")
sp_input.add_argument("--known_dmr", "-kd",
                  action="store",
                  type=str,
                  required= False,
                  default= os.path.join(os.path.dirname(
                                                os.path.realpath(__file__)
                                                    ),
                                             "Imprinted_DMR_List_V1.tsv"),
                  help="The path to the input file for known imprinted DMRs."
                        "File must have the following information in the following column order: "
                        "chromosome\tstart\tend\tMethylatedAlleleOrigin "
                        "where the methylated allele origin must be either "
                        "maternal or paternal (First row must be header). "
                        "By default, we use iDMR list in repo's patmat directory.")
sp_input.add_argument("--whatshap_vcf", "-wv",
                      action="store",
                      type=str,
                      required=False,
                      default=None,
                      help=("Path to the WhatsHap phased vcf file that is produced from "
                            "phasing input vcf file using nanopore reads via WhatsHap. This can be useful "
                            "when the chromosome-scale phased variants are very sparce. "
                            "File must be sorted and indexed using tabix."))
sp_input.add_argument("--whatshap_block", "-wb",
                      action="store",
                      type=str,
                      required=False,
                      default=None,
                      help=("Path to the WhatsHap block file. This file can be"
                            "created using whatshap stats command. File must be" 
                            "converted to a bed format with chromosome\tstart\tend in "
                            "the first three columns (First row must be header). If no block file is given"
                            " then the assumption is that the last part after : sign "
                            "in the 10th column is the phase set (PS) name and blocks will be"
                            " calculated internaly."))
sp_input.add_argument("--black_list", "-bl",
                      action="store",
                      type=str,
                      required=False,
                      default= None,
                      help="List of regions to ignore phased varinats at them."
                      " Three first columns must be chromosome\tstart\tend."
                      " If black list is given the vcf file must be indexed using tabix.")
sp_input.add_argument("--per_read", "-pr",
                      action="store",
                      type=str,
                      required=False,
                      default= None,
                      help="If it is your second try and you have per "
                      "read info file give the path to the per "
                      "read info file. This will be significantly faster."
                      " This is useful when you want to try different thresholds for options,"
                      " different dmr list, black list, include/exclude indels, and include/exclude supp reads.")
sp_input.add_argument("--hapratio", "-hr",
                      action="store",
                      type=float,
                      required=False,
                      default=0.75,
                      help=("0-1 . Minimmum ratio of variants a read must have from a haplotype"
                            " to assign it to that haplotype. Default is 0.75."))
sp_input.add_argument("--min_base_quality", "-mbq",
                      action="store",
                      type=int,
                      required=False,
                      default=7,
                      help=("Only include bases with phred score higher or"
                            " equal than this option. Default is >=7."))
sp_input.add_argument("--mapping_quality", "-mq",
                      action="store",
                      type=int,
                      required=False,
                      default=20,
                      help=("An integer value to specify thereshold for "
                            "filtering reads based om mapping quality. "
                            "Default is >=20"))
sp_input.add_argument("--min_variant", "-mv",
                      action="store",
                      type=int,
                      required=False,
                      default=1,
                      help=("minimum number of phased variants must a read "
                            "have to be phased. Default= 1"))
sp_input.add_argument("--min_read_number", "-mr",
                      action="store",
                      type=int,
                      required=False,
                      default=2,
                      help=("minimum number of reads to support a variant"
                            " to assign to each haplotype. Default= 2"))
sp_input.add_argument("--min_cg", "-mcg",
                      action="store",
                      type=int,
                      required=False,
                      default= 12,
                      help=("Minimmum number of CpGs an iDMR must have to "
                            " consider it for PofO assignment. Default is 12."))
sp_input.add_argument("--meth_difference", "-md",
                      action="store",
                      type=float,
                      required=False,
                      default= 0.35,
                      help=("Methylation difference cutoff for HP1-HP2 or "
                            "HP2-HP1 CpG methylation. Default is 0.35."))
sp_input.add_argument("--cpg_difference", "-cd",
                      action="store",
                      type=float,
                      required=False,
                      default= 0.1,
                      help=("Cut off for the fraction of CpGs between haplotypes "
                            "must be differentially methylated at an iDMR to "
                            "consider it for PofO assignment. Default is 0.1."))
sp_input.add_argument("--methyl_coverage", "-mcov",
                      action="store",
                      type=int,
                      required=False,
                      default=1,
                      help=("Minimmum Coverage at each CpG site when calculating"
                            " methylation frequency. Default is 1."))
sp_input.add_argument("--threads", "-t",
                      action="store",
                      type=int,
                      required=False,
                      default=4,
                      help="Number of parallel processes. Default is 4.")
sp_input.add_argument("--chunk_size", "-cs",
                      action="store",
                      type=int,
                      required=False,
                      default=100,
                      help=("Chunk per process. Default is 100"))
sp_input.add_argument("--include_supplementary", "-is",
                      action="store_true",
                      required=False,
                      help="Also include supplementary reads (Not recommended).")
sp_input.add_argument("--include_indels", "-ind",
                      action="store_true",
                      required=False,
                      help="Also include indels for read phasing to haplotypes.")
sp_input.add_argument('--version', action='version', version='%(prog)s 1.2.0')
args = parser.parse_args()


if __name__ == '__main__':
    main(args)

