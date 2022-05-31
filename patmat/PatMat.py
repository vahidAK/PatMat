#! /usr/bin/env python3
# coding=utf-8

# Copyright (C) 2020  Vahid Akbari

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

def get_snv_info(feed_list,
                  alignment_file,
                  chrom,
                  include_supp):
    read_HP_list= list()
    samfile = pysam.AlignmentFile(alignment_file, 'rb')
    for info in feed_list:
        HP,position,ref,alt= info
        try:
            try:
                sam_pileup= samfile.pileup(chrom,
                                           position,
                                           position+1,
                                           truncate=True)
            except:
                sam_pileup= samfile.pileup(chrom[3:],
                                           position,
                                           position+1,
                                           truncate=True)
        except:#The cordiniate is not found so ignore it
            warnings.warn("{}:{}-{} does not exist in alignment file."
                          " Skipping it.".format(chrom,
                                                   position,
                                                   position+1))
            continue
        for pileupcolumn in sam_pileup:
            pileupcolumn.set_min_base_quality(0)
            if pileupcolumn.pos == position:
                for pileupread in pileupcolumn.pileups:
                    if not include_supp:
                        if pileupread.alignment.is_supplementary:
                            continue
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
                    if pileupread.alignment.is_reverse:
                        strand = "-"
                    else:
                        strand = "+"
                    key_per_read = (chrom,read_start,read_end,
                                    read_id,strand,flag,read_mq)
                    val= [position,phred]
                    read_base= None
                    if len(ref) == 1 and len(alt)==1: # dealing with mismatches
                        read_base= pileupread.alignment.query_sequence[
                                pileupread.query_position]
                        read_base= read_base.upper()
                    elif len(ref) > 1 and len(alt) == 1: #dealing with deletions
                        if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) -1:
                            read_base= pileupread.alignment.query_sequence[
                                    pileupread.query_position]
                            read_base= read_base.upper()
                        elif pileupread.indel == 0:
                            read_base = pileupread.alignment.query_sequence[
                                            pileupread.query_position:pileupread.query_position+len(ref)]
                            read_base = read_base.upper()
                        else:
                            continue
                    elif len(ref) == 1 and len(alt) > 1: #dealing with insertions
                        if pileupread.indel > 0 and pileupread.indel == len(alt) -1:
                            read_base= pileupread.alignment.query_sequence[
                                            pileupread.query_position:pileupread.query_position+len(alt)]
                            read_base= read_base.upper()
                        elif pileupread.indel == 0:
                            read_base= pileupread.alignment.query_sequence[
                                    pileupread.query_position]
                            read_base= read_base.upper()
                        else:
                            continue
                    if read_base is not None:
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
                        elif HP == '1/2' and read_base == ref:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [ref]))])
                        elif HP == '1/2' and read_base == alt:
                            read_HP_list.append([(*key_per_read,"unphased"),
                                            ':'.join(map(str,val + [alt]))])
                            
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


def FrequencyCalculator(file_path):
    """
    Calculating methylation frequency for each phased methylation call
    file.  The output methylation frequency file include fractional
    methylation which can be used for differential methylation analysis
    and detection of differentially methylated regions (DMR)
    """
    dict_mod = defaultdict(int)
    dict_all = defaultdict(int)
    with open(file_path) as input_file:
        next(input_file)  # To exclude header
        for line in input_file:
            line = line.rstrip().split('\t')
            if line[3] == '+':
                key = tuple(line[0:4])
            else:
                key = tuple([line[0],
                                str(int(line[1])+1),
                                str(int(line[2])+1),
                                line[3]])
            dict_all[key] += 1
            if float(line[5]) > 0:
                dict_mod[key] += 1
    return dict_mod, dict_all

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
        dmr_start= int(line[1])
        dmr_end= int(line[2])
        origin= line[3].lower()
        try:
            records = tb_methylcall.query(dmr_chrom, dmr_start, dmr_end)
        except:
            warnings.warn("iDMR {}:{}-{} was ignored because it does not have any reads in "
                          "methylation call file.".format(dmr_chrom, dmr_start, dmr_end))
            records = "NA"
        if records != "NA":
            for record in records:
                key= (record[4],record[3])
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
        if len(common_cg) < min_cg and abs(diff_cg_hp1 - diff_cg_hp2) / len(common_cg) < cpg_difference:# or abs(hp1_freq - hp2_freq) < 0.1:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg_And_cpg_difference" + '\t' +
                           "Ignored:DidNotMeet_min_cg_And_cpg_difference" +'\n')
            continue
        elif len(common_cg) < min_cg:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg" + '\t' +
                           "Ignored:DidNotMeet_min_cg" +'\n')
            continue
        elif abs(diff_cg_hp1 - diff_cg_hp2) / len(common_cg) < cpg_difference:# or abs(hp1_freq - hp2_freq) < 0.1:
            out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_cpg_difference" + '\t' +
                           "Ignored:DidNotMeet_cpg_difference" +'\n')
            continue
        num_cg_to_add_hp1= (diff_cg_hp1 * hp1_freq)/len(common_cg)
        num_cg_to_add_hp2= (diff_cg_hp2 * hp2_freq)/len(common_cg)
        out_meth.write('\t'.join(line)+'\t'+str(len(common_cg))+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                       '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + str(num_cg_to_add_hp1) + '\t' + 
                       str(num_cg_to_add_hp2)+'\n')  
        if diff_cg_hp1 > diff_cg_hp2:# and hp1_freq > hp2_freq:
            if origin == 'maternal':
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += num_cg_to_add_hp1
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += num_cg_to_add_hp2 #0
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += num_cg_to_add_hp1
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += num_cg_to_add_hp2 #0
            if origin == 'paternal':
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += num_cg_to_add_hp2 #0
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += num_cg_to_add_hp1
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += num_cg_to_add_hp2 #0
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += num_cg_to_add_hp1
        elif diff_cg_hp2 > diff_cg_hp1:# and hp2_freq > hp1_freq:
            if origin == 'maternal':
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += num_cg_to_add_hp2
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += num_cg_to_add_hp1 #0
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += num_cg_to_add_hp2
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += num_cg_to_add_hp1 #0
            if origin == 'paternal':
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += num_cg_to_add_hp1 #0
                chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += num_cg_to_add_hp2
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += num_cg_to_add_hp1 #0
                chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += num_cg_to_add_hp2
    dmr_file.close()
    out_meth.close()
    return chrom_hp_origin_count

def vcf2dict(vcf, sites_to_ignore):
    """
    This function converts the input vcf file to haplotype1 and
    haplotype2 dictionaries to be used for read phasing.
    """
    vcf_dict= defaultdict(dict)
    vcf_file = openfile(vcf)
    for line in vcf_file:
        if line.startswith("#"):
            continue
        line_list = line.rstrip().split('\t')
        line_list[9] = line_list[9].replace(".|0","1|0").replace(".|1","0|1").replace("1|.","1|0").replace("0|.","0|1")
        chrom = line_list[0]
        if (line_list[0], line_list[1]) not in sites_to_ignore:
            if (line_list[9].startswith('1|0') or 
                line_list[9].startswith('0|1')):
                vcf_dict[chrom][line_list[1]] = line_list[9].split(':')[0]
    vcf_file.close()
    return vcf_dict


def strand_vcf2dict_phased(vcf_strand, vcf, sites_to_ignore):
    final_dict= defaultdict(set)
    vcf_dict= vcf2dict(vcf_strand, sites_to_ignore)
    with openfile(vcf) as vf:
        for line in vf:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            if (line[0] in vcf_dict and 
                line[1] in vcf_dict[line[0]]):
                final_dict[line[0]].add((vcf_dict[line[0]][line[1]],
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif (line[9].startswith('0/1') or 
                  line[9].startswith('0|1') or
                  line[9].startswith('1|0')):
                final_dict[line[0]].add(("0/1",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif line[9].startswith('1/2') or line[9].startswith('1|2'):
                final_dict[line[0]].add(("1/2",
                                        int(line[1])-1,
                                        line[4].split(',')[0].upper(),
                                        line[4].split(',')[1].upper()))
    return final_dict
                

def vcf2dict_phased(block_file, 
                    vcf_whats, 
                    vcf_strand,
                    hapRatio,
                    minSNV,
                    sites_to_ignore):
    """
    This function converts the input vcf file to haplotype1 and
    haplotype2 dictionaries to be used for read phasing.
    """
    final_dict= defaultdict(set)
    vcf_dict= vcf2dict(vcf_strand, sites_to_ignore)
    tb_whatsvcf= tabix.open(os.path.abspath(vcf_whats))
    with openfile(vcf_whats) as wv:
        for line in wv:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            if line[9].startswith('1/2') or line[9].startswith('1|2'):
                final_dict[line[0]].add(("1/2",
                                        int(line[1])-1,
                                        line[4].split(',')[0].upper(),
                                        line[4].split(',')[1].upper()))
            elif (line[9].startswith('0/1') or 
                  line[9].startswith('0|1') or
                  line[9].startswith('1|0')):
                if (line[0] in vcf_dict and 
                      line[1] in vcf_dict[line[0]]):
                    final_dict[line[0]].add((vcf_dict[line[0]][line[1]],
                                            int(line[1])-1,
                                            line[3].upper(),
                                            line[4].upper()))
                elif line[9].startswith('0/1'):
                    final_dict[line[0]].add(("0/1",
                                            int(line[1])-1,
                                            line[3].upper(),
                                            line[4].upper()))
    for blocks in block_file:
        b_chrom, b_start, b_end= blocks
        agreement_count= 0
        disagreement_count= 0
        try:
            records_whats = tb_whatsvcf.query(b_chrom, b_start-1, b_end+1)
            records_whats = list(records_whats)
        except:
            warnings.warn("{}:{}-{} block does not exist in the "
                          "vcf file. Skipping it.".format(b_chrom, b_start-1, b_end+1))
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
                pos = int(vcf_line[1])-1#VCF file is 1-based
                if (vcf_line[9].startswith('1|0') or 
                    vcf_line[9].startswith('0|1')):
                    if (agreement_count > disagreement_count and
                        agreement_count >= minSNV and 
                        agreement_count/(agreement_count+disagreement_count) >= hapRatio and
                        vcf_line[1] not in disagreement_sites):
                        final_dict[b_chrom].add((vcf_line[9].split(':')[0],
                                                pos,
                                                vcf_line[3].upper(),
                                                vcf_line[4].upper()))
                    elif (disagreement_count > agreement_count and
                          disagreement_count >= minSNV and 
                          disagreement_count/(agreement_count+disagreement_count) >= hapRatio and
                          vcf_line[1] not in agreement_sites):
                        final_dict[b_chrom].add((vcf_line[9].split(':')[0][::-1],
                                                pos,
                                                vcf_line[3].upper(),
                                                vcf_line[4].upper()))
                    else:
                        if vcf_line[1] not in vcf_dict[vcf_line[0]]: 
                            final_dict[b_chrom].add(("0/1",
                                                    pos,
                                                    vcf_line[3].upper(),
                                                    vcf_line[4].upper()))
        else:
            for vcf_line in records_whats:
                pos = int(vcf_line[1])-1#VCF file is 1-based
                if (vcf_line[9].startswith('1|0') or 
                    vcf_line[9].startswith('0|1')):
                    final_dict[b_chrom].add(("0/1",
                                            pos,
                                            vcf_line[3].upper(),
                                            vcf_line[4].upper()))
    vcf_dict.clear()
    return final_dict


def per_read_snv(vcf_dict,
                 bam_file,
                 chunk,
                 threads,
                 include_supplementary,
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
            description= "Tagging SNVs to reads from {}: ".format(chrom)
            with tqdm(total=len(feed_list),
                desc=description,
                bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                                  ) as pbar:
                for vcf_info_list in feed_list:
                    p= mp.Pool(len(vcf_info_list))
                    results= p.starmap(get_snv_info,
                                       list(zip(vcf_info_list,
                                                repeat(bam_file),
                                                repeat(chrom),
                                                repeat(include_supplementary))))
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
                    hp1_snvs= per_read_hp[key]['HP1']
                    hp1_count= str(len(hp1_snvs))
                except:
                    hp1_count= '0'
                if hp1_count == '0':
                    hp1_snvs= ['NA']
                    
                try:
                    hp2_snvs= per_read_hp[key]['HP2']
                    hp2_count= str(len(hp2_snvs))
                except:
                    hp2_count= '0'   
                if hp2_count == '0':
                    hp2_snvs= ['NA']
                    
                try:
                    unphased_snvs= per_read_hp[key]['unphased']
                    unphased_count= str(len(unphased_snvs))
                except:
                    unphased_count= '0'
                if unphased_count == '0':
                    unphased_snvs= ['NA']
                    
                out_to_write= list(map(str,key)) + [','.join(hp1_snvs), 
                                                    ','.join(hp2_snvs),
                                                    ','.join(unphased_snvs)]
                
                perReadinfo.write('\t'.join(out_to_write)+'\n')
        else:
            warnings.warn("{} does not have any mapped reads in alignment "
                          "file Or alignment is truncated or corrupt indexed. "
                          "Skipping it.".format(chrom))
            # unprocessed_chrom.add(chrom)
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



def main_phase(args):
    """
    This is the phase module which phase the nanopore reads and
    methylation data to corresponding haplotype using vcf file and
    processed methylation call file.
    """
    hapRatio = args.hapratio
    minSNV= args.min_snv
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
                records= tb_vcf.query(line[0], int(line[1]), int(line[2])+1)
                for record in records:
                    sites_to_ignore.add((record[0],record[1]))
    if args.known_dmr is not None:
        try:
            known_dmr= os.path.abspath(args.known_dmr)
            MethylCallfile = os.path.abspath(args.methylcallfile)
            tb_methylcall = tabix.open(MethylCallfile)
        except:
            raise Exception("When PofO needs to be determined, "
                            "known DMRs and NanoMethPhased processed methylationCall file"
                            " must be specified")
        if not os.path.isfile(MethylCallfile+".tbi"):
            raise Exception("It seems processed methylation call file "
                            "is not index by tabix.")
    if args.per_read is not None:
        per_read_file= os.path.abspath(args.per_read)
    else:
        if args.strand_vcf is not None and args.whatshap_vcf is None:
            warnings.warn("No WhatsHap phased vcf is given. Using strand-seq phased vcf only.")
            vcf_strand = os.path.abspath(args.strand_vcf)
            final_dict= strand_vcf2dict_phased(vcf_strand, vcf, sites_to_ignore)
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
            warnings.warn("Using both strand-seq phased and WhatsHap phased vcf.")
            vcf_whats= os.path.abspath(args.whatshap_vcf)
            if not os.path.isfile(MethylCallfile+".tbi"):
                raise Exception("It seems that whatshap vcf "
                                "is not index by tabix.")
            vcf_strand = os.path.abspath(args.strand_vcf)
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
                                        minSNV,
                                        sites_to_ignore)
        else:
            raise Exception("No strand-seq and/or WhatsHap phased vcf is given.")
    if args.per_read is None:
        out_per_read = out + '_HP1'
        perReadinfo= open(out_per_read+"_HP2_PerReadInfo.tsv", 'w')
        perReadinfo.write("#Chromosome\tReadRefStart\tReadRefEnd\tReadID\t"
                          "Strand\tReadFlag\tReadMapQuality\t"
                          "Position:BaseQuality:AltAllele_HP1-SNVs\t"
                          "Position:BaseQuality:AltAllele_HP2-SNVs\t"
                          "Position:BaseQuality:AltAllele_UnPhasedHetSNVs\n")
        per_read_snv(final_dict,
                     bam_file,
                     chunk,
                     threads,
                     args.include_supplementary,
                     perReadinfo)
        final_dict.clear()
        perReadinfo.close()
        per_read= openfile(out_per_read+"_HP2_PerReadInfo.tsv")
    else:
        per_read= openfile(per_read_file)
    chrom_list= set()
    read_Unphased = defaultdict(set)
    read_dictHP1 = defaultdict(set)
    read_dictHP2 = defaultdict(set)
    all_reads= defaultdict(set)
    for line in per_read:
        if line.startswith("#"):
            continue
        line= line.rstrip().split('\t')
        if int(line[6]) < MappingQuality:
            continue
        chrom_list.add(line[0])
        key= tuple(line[3:5])
        hp1s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[7].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality]
        hp2s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[8].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality]
        unassigns= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[9].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality]
        all_reads[line[0]].add(key)
        hp1_count= len(hp1s)
        hp2_count= len(hp2s)
        if (hp1_count > hp2_count and 
            hp1_count/(hp1_count+hp2_count) >= hapRatio and 
            hp1_count >= minSNV):
            read_dictHP1[line[0]].add(key)
            
        elif (hp2_count > hp1_count and 
              hp2_count/(hp1_count+hp2_count) >= hapRatio and 
              hp2_count >= minSNV):
            read_dictHP2[line[0]].add(key)
        else:
            read_Unphased[line[0]].add(key)
    per_read.close()
    all_reads.clear()
    snv_dict_HP1= defaultdict(lambda: defaultdict(int))
    snv_dict_HP2= defaultdict(lambda: defaultdict(int))
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
        key= tuple(line[3:5])
        hp1s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[7].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality]
        hp2s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[8].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality]
        unassigns= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[9].split(',') if i != 'NA' and 
                      int(i.split(':')[1]) >= MinBaseQuality]
        if key in read_dictHP1[line[0]]:
            for i in hp1s+hp2s+unassigns:
                snv_dict_HP1[(i[0],i[1])][i[2]] += 1
                snv_dict_HP2[(i[0],i[1])][i[2]] += 0
        elif key in read_dictHP2[line[0]]:
            for i in hp1s+hp2s+unassigns:
                snv_dict_HP1[(i[0],i[1])][i[2]] += 0
                snv_dict_HP2[(i[0],i[1])][i[2]] += 1
    per_read.close()
    # if args.per_read is not None:
    #     per_read= openfile(os.path.abspath(args.per_read))
    # else:
    #     per_read= openfile(out_per_read+"_HP2_PerReadInfo.tsv")
    # for line in per_read:
    #     if line.startswith("#"):
    #         continue
    #     line= line.rstrip().split('\t')
    #     if int(line[6]) < MappingQuality:
    #         continue
    #     key= tuple(line[3:5])
    #     if key in read_Unphased[line[0]]:
    #         hp1s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[7].split(',') if i != 'NA' and 
    #                       int(i.split(':')[1]) >= MinBaseQuality]
    #         hp2s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[8].split(',') if i != 'NA' and 
    #                       int(i.split(':')[1]) >= MinBaseQuality]
    #         unassigns= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[9].split(',') if i != 'NA' and 
    #                       int(i.split(':')[1]) >= MinBaseQuality]
    #         for i in hp1s+hp2s+unassigns:
    #             if (i[0],i[1]) in snv_dict_HP1 and i[2] in snv_dict_HP1[(i[0],i[1])]:
    #                 hp1_count = snv_dict_HP1[(i[0],i[1])][i[2]]
    #             if (i[0],i[1]) in snv_dict_HP2 and i[2] in snv_dict_HP2[(i[0],i[1])]:
    #                 hp2_count = snv_dict_HP2[(i[0],i[1])][i[2]]
    #         if (hp1_count > hp2_count and 
    #             hp1_count/(hp1_count+hp2_count) >= hapRatio and 
    #             hp1_count >= minSNV):
    #             read_dictHP1[line[0]].add(key)
    #         elif (hp2_count > hp1_count and 
    #               hp2_count/(hp1_count+hp2_count) >= hapRatio and 
    #               hp2_count >= minSNV):
    #             read_dictHP2[line[0]].add(key)
    # per_read.close()
    if not read_dictHP1 and not read_dictHP2:
        raise Exception("No reads could be phased. Probably phased vcf files"
                        " do not have any phased SNVs or have very few.")
    out_non_pofo = out + '_NonPofO_Hp1-HP2_Reassignment.vcf'
    out_non_pofo_reads = out + '_NonPofO_Hp1-HP2_reads.tsv'
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
            if line[0] in chrom_list:
                if (line[9].startswith("0/1") or
                    line[9].startswith("0|1") or
                    line[9].startswith("1|0")):
                    hp1_count_alt= snv_dict_HP1[(line[0],str(int(line[1])-1))][line[4]]
                    hp2_count_alt= snv_dict_HP2[(line[0],str(int(line[1])-1))][line[4]]
                    hp1_count_ref= snv_dict_HP1[(line[0],str(int(line[1])-1))][line[3]]
                    hp2_count_ref= snv_dict_HP2[(line[0],str(int(line[1])-1))][line[3]]
                    if ((hp1_count_alt > hp2_count_alt and
                         hp1_count_alt > hp1_count_ref and
                         hp1_count_alt >= min_read_reassignment) or
                        (hp2_count_ref > hp1_count_ref and
                         hp2_count_ref > hp2_count_alt and
                         hp2_count_ref >= min_read_reassignment)):
                        hp1s.add((line[0],line[1],line[4]))
                        hp2s.add((line[0],line[1],line[3]))
                        re_assignment.write('\t'.join(line[0:-2]+
                                                      [':'.join(line[-2].split(":")).replace(":PS", "")+":PS"]+
                                                      ["0|1:"+':'.join(line[-1].split(":")[1:])+":HP2|HP1"])+'\n')
                    elif ((hp2_count_alt > hp1_count_alt and
                           hp2_count_alt > hp2_count_ref and
                           hp2_count_alt >= min_read_reassignment) or
                          (hp1_count_ref > hp2_count_ref and
                           hp1_count_ref > hp1_count_alt and
                           hp1_count_ref >= min_read_reassignment)):
                        hp1s.add((line[0],line[1],line[3]))
                        hp2s.add((line[0],line[1],line[4]))
                        re_assignment.write('\t'.join(line[0:-2]+
                                                      [':'.join(line[-2].split(":")).replace(":PS", "")+":PS"]+
                                                      ["1|0:"+':'.join(line[-1].split(":")[1:])+":HP1|HP2"])+'\n')
                    else:
                        re_assignment.write('\t'.join(line[0:-2]+
                                                      [':'.join(line[-2].split(":")).replace(":PS", "")]+
                                                      [':'.join(line[-1].split(":")).replace("|", "/")])+'\n')
                elif line[9].startswith("1/1") or line[9].startswith("1|1"):
                    re_assignment.write('\t'.join(line)+'\n')
                
                elif line[9].startswith("1/2") or line[9].startswith("1|2"):
                    hp1_count_alt= snv_dict_HP1[(line[0],str(int(line[1])-1))][line[4].split(',')[1].upper()]
                    hp2_count_alt= snv_dict_HP2[(line[0],str(int(line[1])-1))][line[4].split(',')[1].upper()]
                    hp1_count_ref= snv_dict_HP1[(line[0],str(int(line[1])-1))][line[4].split(',')[0].upper()]
                    hp2_count_ref= snv_dict_HP2[(line[0],str(int(line[1])-1))][line[4].split(',')[0].upper()]
                    if ((hp1_count_alt > hp2_count_alt and
                         hp1_count_alt > hp1_count_ref and
                         hp1_count_alt >= min_read_reassignment) or
                        (hp2_count_ref > hp1_count_ref and
                         hp2_count_ref > hp2_count_alt and
                         hp2_count_ref >= min_read_reassignment)):
                        hp1s.add((line[0],line[1],line[4].split(',')[1]))
                        hp2s.add((line[0],line[1],line[4].split(',')[0]))
                        re_assignment.write('\t'.join(line[0:4]+[line[4].split(',')[0]+","+line[4].split(',')[1]]+
                                                      line[5:-2]+
                                                      [':'.join(line[-2].split(":")).replace(":PS", "")+":PS"]+
                                                      ["1|2:"+':'.join(line[-1].split(":")[1:])+":Ref/HP2|HP1"])+'\n')
                    elif ((hp2_count_alt > hp1_count_alt and
                           hp2_count_alt > hp2_count_ref and
                           hp2_count_alt >= min_read_reassignment) or
                          (hp1_count_ref > hp2_count_ref and
                           hp1_count_ref > hp1_count_alt and
                           hp1_count_ref >= min_read_reassignment)):
                        hp1s.add((line[0],line[1],line[4].split(',')[0]))
                        hp2s.add((line[0],line[1],line[4].split(',')[1]))
                        re_assignment.write('\t'.join(line[0:4]+[line[4].split(',')[0]+","+line[4].split(',')[1]]+
                                                      line[5:-2]+
                                                      [':'.join(line[-2].split(":")).replace(":PS", "")+":PS"]+
                                                      ["1|2:"+':'.join(line[-1].split(":")[1:])+":Ref/HP1|HP2"])+'\n')
                    else:
                        re_assignment.write('\t'.join(line[0:-2]+
                                                      [':'.join(line[-2].split(":")).replace(":PS", "")]+
                                                      [':'.join(line[-1].split(":")).replace("|", "/")])+'\n')
                else:
                    re_assignment.write('\t'.join(line)+'\n')
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
        key= tuple(line[3:5])
        variants= [(line[0],str(int(i.split(":")[0])+1),i.split(':')[2]) for i in 
                    set(line[7].split(',') + line[8].split(',') +line[9].split(',')) 
                    if i != 'NA' and int(i.split(':')[1]) >= MinBaseQuality]
        hp1_count= len([i for i in variants if i in hp1s])
        hp2_count= len([i for i in variants if i in hp2s])
        if (hp1_count > hp2_count and 
            hp1_count/(hp1_count+hp2_count) >= hapRatio and 
            hp1_count >= minSNV):
            read_dictHP1[line[0]].add(key)
        elif (hp2_count > hp1_count and 
              hp2_count/(hp1_count+hp2_count) >= hapRatio and 
              hp2_count >= minSNV):
            read_dictHP2[line[0]].add(key)
    per_read.close()
    for chrom, reads in read_dictHP1.items():
        for read in reads:
            reads_NonPofO.write(chrom+'\t'+'\t'.join(read)+'\t'+"HP1"+'\n')
    for chrom, reads in read_dictHP2.items():
        for read in reads:
            reads_NonPofO.write(chrom+'\t'+'\t'.join(read)+'\t'+"HP2"+'\n')
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
                        reads_PofO.write(chrom+'\t'+'\t'.join(read)+'\t'+"maternal"+'\n')
                if chrom_hp_origin[chrom]['HP1'][0] == 'paternal':
                    for read in reads:
                        reads_PofO.write(chrom+'\t'+'\t'.join(read)+'\t'+"paternal"+'\n')
            for chrom, reads in read_dictHP2.items():
                if chrom not in chrom_hp_origin:
                    continue
                if chrom_hp_origin[chrom]['HP2'][0] == 'maternal':
                    for read in reads:
                        reads_PofO.write(chrom+'\t'+'\t'.join(read)+'\t'+"maternal"+'\n')
                if chrom_hp_origin[chrom]['HP2'][0] == 'paternal':
                    for read in reads:
                        reads_PofO.write(chrom+'\t'+'\t'.join(read)+'\t'+"paternal"+'\n')
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
                                assignment_file.write('\t'.join(line[0:-1])+'\t'+
                                                          line[-1].replace("HP1","Mat").replace("HP2","Pat")
                                                          +'\n')
                            if chrom_hp_origin[line[0]]['HP1'][0] == 'paternal':
                                assignment_file.write('\t'.join(line[0:-1])+'\t'+
                                                          line[-1].replace("0|1","1|0").replace("HP1","Pat").replace("HP2","Mat")
                                                          +'\n')
                        elif line[9].startswith('1|0'):
                            if chrom_hp_origin[line[0]]['HP2'][0] == 'maternal':
                                assignment_file.write('\t'.join(line[0:-1])+'\t'+
                                                          line[-1].replace("1|0","0|1").replace("HP1","Pat").replace("HP2","Mat")
                                                          +'\n')
                            if chrom_hp_origin[line[0]]['HP2'][0] == 'paternal':
                                assignment_file.write('\t'.join(line[0:-1])+'\t'+
                                                          line[-1].replace("HP1","Mat").replace("HP2","Pat")
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
                            assignment_file.write('\t'.join(line[0:-2]+
                                                            [':'.join(line[-2].split(":")).replace(":PS", "")]+
                                                            [':'.join(line[-1].split(":")).replace("|", "/").
                                                             replace(":Ref/HP1/HP2","").
                                                             replace(":HP1/HP2","").
                                                             replace(":HP2/HP1","")])+'\n')
                    else:
                        assignment_file.write('\t'.join(line[0:-2]+
                                                        [':'.join(line[-2].split(":")).replace(":PS", "")]+
                                                        [':'.join(line[-1].split(":")).replace("|", "/").
                                                          replace(":Ref/HP1/HP2","").
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
 
    
def phase_parser(subparsers):
    """
    Specific argument parser for phase command.
    """
    sub_phase = subparsers.add_parser("phase",
                                      add_help=False,
                                      help="Phasing reads and Methylation.",
                                      description="Phasing reads and "
                                      "Methylation")
    sp_input = sub_phase.add_argument_group("required arguments")
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
                          help="The path to the chromosome-scale Strand-seq phased vcf file. ")
                           # "If it is your second try and you have per read "
                           # "info file from the first try there is no need to "
                           # "give vcf file, instead give the path to the per "
                           # "read info file using --per_read option which will "
                           # "be significantly faster.")
    sp_input = sub_phase.add_argument_group("required arguments if PofO needs to be determined.")
    sp_input.add_argument("--known_dmr", "-kd",
                      action="store",
                      type=str,
                      required= False,
                      default= os.path.join(os.path.dirname(
                                                    os.path.realpath(__file__)
                                                        ),
                                                 "Imprinted_DMR_List_V1.tsv"),
                      help="The path to the input file for known imprinted DMRs."
                      "File must have the following information the following column order: "
                      "chromosome\tstart\tend\tMethylatedAlleleOrigin "
                      "where origine is the methylated allele origine which must be either "
                      "maternal or paternal. By default, we use version 1 list in repo's patmat directory.")
    sp_input.add_argument("--methylcallfile", "-mc",
                          action="store",
                          type=str,
                          required=False,
                          default=None,
                          help=("If you want to phase methyl call file "
                                "(methycall output format) to also calculate "
                                "methylation frequency for each haplotype "
                                "give the path to the bgziped methylation "
                                "call file from methyl_call_processor Module."
                                ))
    sp_input = sub_phase.add_argument_group("Optional arguments.")
    sp_input.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    sp_input.add_argument("--whatshap_vcf", "-wv",
                          action="store",
                          type=str,
                          required=False,
                          default=None,
                          help=("Path to the WhatsHap phased vcf file that is produced from "
                                "phasing nanopore reads using WhatsHap. This can be useful "
                                "when the chromosome-scale phased variants are very sparce. "
                                "File must be sorted and indexed using tabix."))
    sp_input.add_argument("--whatshap_block", "-wb",
                          action="store",
                          type=str,
                          required=False,
                          default=None,
                          help=("Path to the WhatsHap block file file. This file can be"
                                "created using whatshap stats command. File must be" 
                                "converted to a bed format with chromosome\tstart\tend in "
                                "the first three columns. If no block file is given"
                                " then the assumption is that the last part after : sign "
                                "in the 10th column is the phase set (PS) name and blocks will be"
                                " calculated internaly."))
    sp_input.add_argument("--black_list", "-bl",
                          action="store",
                          type=str,
                          required=False,
                          default= None,
                          help="List of regions to ignore ther strand-seq phasedcstatus"
                          " three first columns must be chromosome\tstart\tend."
                          " If black list is given the vcf file must be indexed using tabix.")
    sp_input.add_argument("--per_read", "-pr",
                          action="store",
                          type=str,
                          required=False,
                          default= None,
                          help="If it is your second try and you have per "
                          "read info file give the path to the per "
                          "read info file. This will be significantly faster.")
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
    sp_input.add_argument("--min_snv", "-ms",
                          action="store",
                          type=int,
                          required=False,
                          default=1,
                          help=("minimum number of phased SNVs must a read "
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
                          default= 11,
                          help=("Minimmum number of CpGs an iDMR must have to "
                                " consider it for PofO assignment. Default is 11."))
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
                          help="Also include supplementary reads.")
    sub_phase.set_defaults(func=main_phase)


def main():
    """
    Docstring placeholder.
    """
    parser = argparse.ArgumentParser(
        prog="PatMat.py",
        description="PatMat")
    subparsers = parser.add_subparsers(title="Modules")
    phase_parser(subparsers)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()

