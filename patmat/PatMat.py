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

import subprocess
import os
import gzip
import bz2
import argparse
import warnings
import multiprocessing as mp
from collections import defaultdict
from itertools import repeat
import pysam
import tabix
import re
import math
from tqdm import tqdm
from modbampy import ModBam

def get_variant_info(feed_list,
                  alignment_file,
                  chrom):
    read_HP_list= list()
    samfile = pysam.AlignmentFile(alignment_file, 'rb')
    NA_phred= False
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
                    read_len= pileupread.alignment.query_alignment_length
                    read_start= pileupread.alignment.reference_start
                    read_end= pileupread.alignment.reference_end
                    flag= pileupread.alignment.flag
                    try:
                        phred= pileupread.alignment.query_qualities[
                                                    pileupread.query_position]
                    except:
                        phred= "NA"
                        NA_phred= True
                        
                    if pileupread.alignment.is_supplementary:
                        suppl_flag= str(flag)+":yes"
                    else:
                        suppl_flag= str(flag)+":no"
                    if pileupread.alignment.is_reverse:
                        strand = "-"
                    else:
                        strand = "+"
                    key_per_read = (chrom,read_start,read_end,
                                    read_id,strand,suppl_flag,
                                    str(read_len)+":"+str(read_mq))
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
    return read_HP_list, NA_phred



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
                 out,
                 min_cg,
                 cpg_difference):
    tb_calldml= tabix.open(out+"_callDML.tsv.gz")
    tb_dmltest= tabix.open(out+"_DMLtest.tsv.gz")
    out_meth= open(out+"_CpG-Methylation-Status-at-DMRs.tsv",'w')
    dmr_file= openfile(known_dmr)
    header= next(dmr_file).rstrip()
    out_meth.write(header+"\tAll_CpGs_At_iDMR_CouldBeExaminedInBothHaplotypes\tDifferentiallyMethylatedCpGs_HypermethylatedOnHP1\t"
                    "DifferentiallyMethylatedCpGs_HypermethylatedOnHP2\tMethylationFrequency_HP1\tMethylationFrequency_HP2\t"
                    "DetectionValue_HP1\tDetectionValue_HP2\n")
    chrom_hp_origin_count= defaultdict(lambda: defaultdict(int))
    for line in dmr_file:
        line= line.rstrip().split('\t')
        dmr_chrom= line[0]
        if dmr_chrom not in chrom_list:
            continue
        dmr_start= int(line[1]) -1
        dmr_end= int(line[2])
        origin= line[3].lower()
        if origin not in ["maternal","paternal"]:
            warnings.warn("iDMR: {} does not have a valid origin (must be paternal or maternal. case insensitive). "
                          "This iDMR will not be used for PofO assignment and PofO"
                          " assignment score calculation.".format('\t'.join(line)))
        try:
            records_all = tb_dmltest.query(dmr_chrom, dmr_start, dmr_end)
        except:
            warnings.warn("iDMR {}:{}-{} was ignored because it does not have any CpG in "
                          "DMLtest file.".format(dmr_chrom, dmr_start, dmr_end))
            records_all = "NA"
        try:
            records = tb_calldml.query(dmr_chrom, dmr_start, dmr_end)
        except:
            warnings.warn("iDMR {}:{}-{} was ignored because it does not have any CpG in "
                          "DMLtest file.".format(dmr_chrom, dmr_start, dmr_end))
            records = "NA"
        num_cg = 0
        hp1_freq= 0
        hp2_freq= 0
        diff_cg_hp1 = 0
        diff_cg_hp2 = 0
        if records_all != "NA":
            for record_all in records_all:
                num_cg += 1
                hp1_freq += float(record_all[3])
                hp2_freq += float(record_all[4])
        if records != "NA":
            for record in records:
                if (float(record[3]) - float(record[4])) > 0:
                    diff_cg_hp1 += 1
                elif (float(record[4]) - float(record[3])) > 0:
                    diff_cg_hp2 += 1
        if num_cg < 1:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+"NA"+'\t'+"NA"+
                            '\t'+"NA" + '\t' + "NA"+'\t' + "Ignored:No_CpG" + '\t' +
                            "Ignored:No_CpG" +'\n')
            continue
        hp1_freq= hp1_freq/num_cg
        hp2_freq= hp2_freq/num_cg
        if num_cg < min_cg and abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                            '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg_And_cpg_difference" + '\t' +
                            "Ignored:DidNotMeet_min_cg_And_cpg_difference" +'\n')
            continue
        elif num_cg < min_cg and diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                            '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" + '\t' +
                            "Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" +'\n')
            continue
        elif num_cg < min_cg:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                            '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_min_cg" + '\t' +
                            "Ignored:DidNotMeet_min_cg" +'\n')
            continue
        elif abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                            '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DidNotMeet_cpg_difference" + '\t' +
                            "Ignored:DidNotMeet_cpg_difference" +'\n')
            continue
        elif diff_cg_hp1 == 0 and diff_cg_hp2 == 0: 
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+ 
                            '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + "Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" + '\t' +
                            "Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes" +'\n') 
            continue
        num_cg_to_add_hp1= diff_cg_hp1
        num_cg_to_add_hp2= diff_cg_hp2
        out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
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


def out_freq(chromosome,
             methfile,
             tool,
             hp1_dict,
             hp2_dict,
             callthresh, 
             reference):
    hp1_freq= defaultdict(list)
    hp2_freq= defaultdict(list)
    if tool=="nanopolish":
        chrom_index= 0
        readID_index= 4
        start_index= 2
        strand_index= 1
        modprob_index= 5
    elif tool=="deepsignal":
        chrom_index= 0
        readID_index= 4
        start_index= 1
        strand_index= 2
        modprob_index= 7
    elif tool=="megalodon":
        chrom_index= 1
        readID_index= 0
        start_index= 3
        strand_index= 2
        modprob_index= 4
    if tool=="megalodon" or tool=="deepsignal":     
        with openfile(methfile) as methcall:
            next(methcall)
            for line in methcall:
                line=line.rstrip().split('\t')
                chrom= line[chrom_index]
                if chromosome != chrom:
                    continue
                read_id= line[readID_index]
                strand= line[strand_index]
                cpg_pos= int(line[start_index])
                if strand == "-":
                    cpg_pos= cpg_pos - 1
                if tool=="megalodon":
                    deltaprob= math.exp(float(line[4])) - (1 - math.exp(float(line[4])))
                elif tool=="deepsignal":
                    deltaprob = float(line[7]) - float(line[6])
                if abs(deltaprob) < callthresh:
                    continue
                if chrom in hp1_dict and (chrom,read_id,strand) in hp1_dict[chrom]:
                    if deltaprob > 0:
                        hp1_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(1)
                    else:
                        hp1_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(0)
                if chrom in hp2_dict and (chrom,read_id,strand) in hp2_dict[chrom]:
                    if deltaprob > 0:
                        hp2_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(1)
                    else:
                        hp2_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(0)
    elif tool=="nanopolish":     
        with openfile(methfile) as methcall:
            next(methcall)
            for line in methcall:
                line=line.rstrip().split('\t')
                chrom= line[chrom_index]
                if chromosome != chrom:
                    continue
                read_id= line[readID_index]
                strand= line[strand_index]
                cpg_pos= int(line[start_index])
                logratio = float(line[5])/int(line[9])
                sequence= line[10].upper()
                if abs(logratio) < callthresh:
                    continue
                if chrom in hp1_dict and (chrom,read_id,strand) in hp1_dict[chrom]:
                    if logratio > 0:
                        hp1_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(1)
                    else:
                        hp1_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(0)
                    if int(line[9]) > 1:  # Check if the line includes multi-group CpG
                        splited_groupIndexes = [(j.start())
                                                for j in re.finditer("CG", sequence)]
                        for splited_groupIndex in splited_groupIndexes[1:]:
                            position = cpg_pos + (splited_groupIndex
                                                     - splited_groupIndexes[0])
                            if logratio > 0:
                                hp1_freq[(chrom,str(position),str(position+1))].append(1)
                            else:
                                hp1_freq[(chrom,str(position),str(position+1))].append(0)
                if chrom in hp2_dict and (chrom,read_id,strand) in hp2_dict[chrom]:
                    if logratio > 0:
                        hp2_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(1)
                    else:
                        hp2_freq[(chrom,str(cpg_pos),str(cpg_pos+1))].append(0)
                    if int(line[9]) > 1:  # Check if the line includes multi-group CpG
                        splited_groupIndexes = [(j.start())
                                                for j in re.finditer("CG", sequence)]
                        for splited_groupIndex in splited_groupIndexes[1:]:
                            position = cpg_pos + (splited_groupIndex
                                                     - splited_groupIndexes[0])
                            if logratio > 0:
                                hp2_freq[(chrom,str(position),str(position+1))].append(1)
                            else:
                                hp2_freq[(chrom,str(position),str(position+1))].append(0)
    elif tool=="guppy":
        fasta= pysam.Fastafile(reference)
        chrom_seq= fasta.fetch(reference=chromosome)
        fasta.close()
        with ModBam(methfile) as bam:
            try:
                bamiter= bam.reads(chromosome, 0, len(chrom_seq))
            except:
                warnings.warn("Chromosome {} from reference does no exist in "
                              "guppy bam file. Skipping it.".format(chromosome))
                bamiter= None
            if bamiter is not None:
                for read in bamiter:
                    for pos_mod in read.mod_sites:
                        base_pos= pos_mod.rpos
                        strand= pos_mod.strand
                        if strand=="-":
                            base_pos = base_pos-1
                        if base_pos < 0 or pos_mod.mbase != "m":
                            continue
                        if chrom_seq[base_pos].upper() != 'C' or chrom_seq[base_pos+1].upper() != 'G':
                            continue
                        read_id= pos_mod.query_name
                        cpg_pos= base_pos
                        deltaprob= (pos_mod.qual/256) - (1-(pos_mod.qual/256))
                        if abs(deltaprob) < callthresh:
                            continue
                        if (chromosome,read_id,strand) in hp1_dict[chromosome]:
                            if deltaprob > 0:
                                hp1_freq[(chromosome,str(cpg_pos),str(cpg_pos+1))].append(1)
                            else:
                                hp1_freq[(chromosome,str(cpg_pos),str(cpg_pos+1))].append(0)
                        if (chromosome,read_id,strand) in hp2_dict[chromosome]:
                            if deltaprob > 0:
                                hp2_freq[(chromosome,str(cpg_pos),str(cpg_pos+1))].append(1)
                            else:
                                hp2_freq[(chromosome,str(cpg_pos),str(cpg_pos+1))].append(0)
    else:
        raise Exception("Select tool and call threshold currectly.")
    return (hp1_freq, hp2_freq)
                

def out_pofo_freq(hp_fre,
                  chrom,
                  output):
    if os.path.exists(hp_fre):
        with openfile(hp_fre) as file:
            next(file)
            for line in file:
                line=line.rstrip().split('\t')
                if line[0]==chrom:
                    output.write('\t'.join(line)+'\n')
    else:
        with openfile(hp_fre+".gz") as file:
            next(file)
            for line in file:
                line=line.rstrip().split('\t')
                if line[0]==chrom:
                    output.write('\t'.join(line)+'\n')


def vcf2dict(vcf_strand):
    """
    This function converts the input vcf file to haplotype1 and
    haplotype2 dictionaries to be used for read phasing.
    """
    strand_vars= 0
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
            strand_vars += 1
            vcf_dict[chrom][line_list[1]] = line_list[9].split(':')[0]
    vcf_file.close()
    return vcf_dict, strand_vars


def check_vcfs(vcf,
           vcf_strand,
           vcf_whats):
    common_strand= 0
    common_whats= 0
    vcf_dict, strand_vars= vcf2dict(vcf_strand)
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
                common_strand += 1
    if vcf_whats is None:
        warnings.warn("Out of {} reference heterozygous phased variants in strand-seq vcf, "
                      "{} ({})% variants had the same position in input "
                      "vcf file. These numbers help you to check if"
                      " the input vcf and strand-seq vcf belong to the same sample."
                      "".format(strand_vars, 
                                common_strand,
                                round((common_strand/strand_vars)*100,2)))
    if vcf_whats is not None:
        with openfile(vcf_whats) as vf_w:
            for line in vf_w:
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
                    common_whats += 1
        warnings.warn("Out of {} reference heterozygous phased variants in strand-seq vcf,"
                      " {} ({})% variants had the same position in "
                      "input vcf and {} ({})% variants "
                      " had the same position in whatshap vcf file. "
                      "These numbers help you to check if the input "
                      "vcf, whatshap vcf and strand-seq vcf belong to the same sample."
                      "".format(strand_vars, 
                                common_strand,
                                round((common_strand/strand_vars)*100,2),
                                common_whats,
                                round((common_whats/strand_vars)*100,2)))
    
        
def strand_vcf2dict_phased(vcf_strand,
                           vcf):
    check_vcfs(vcf, vcf_strand, None)
    final_dict= defaultdict(set)
    vcf_dict, strand_vars= vcf2dict(vcf_strand)
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
    check_vcfs(vcf, vcf_strand, vcf_whats)
    final_dict= defaultdict(set)
    vcf_dict, strand_vars= vcf2dict(vcf_strand)
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
                phred_check= False
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
                            if result[1]:
                                phred_check= True
                            for read_info in result[0]:
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
            if phred_check:
                warnings.warn("Some or all bases in some or all reads from {} do "
                              "not have based qualities in the bam file."
                              " Phred quality thereshold will be ignored for these bases.".format(chrom))
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
            # Non-reference het indels (e.g. 1/2) are ignored as they will be in unphased column of per-read info file
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
    cpg_difference= args.cpg_difference
    MinBaseQuality = args.min_base_quality
    chunk = args.chunk_size
    min_read_reassignment= args.min_read_number
    out = os.path.abspath(args.output)
    MappingQuality = args.mapping_quality
    vcf= os.path.abspath(args.vcf)
    known_dmr= os.path.abspath(args.known_dmr)
    MethylCallfile = os.path.abspath(args.methylcallfile)
    ref_file= os.path.abspath(args.reference)
    tool,callthresh= args.tool_and_callthresh.split(':')
    tool= tool.lower()
    callthresh= float(callthresh)
    if tool not in ["guppy","nanopolish","megalodon","deepsignal"]:
        raise Exception("Select tool_and_callthresh option correctly.")
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
        if int(line[6].split(':')[1]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        chrom_list.add(line[0])
        key= (line[0],line[3],line[4])
        hp1s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[7].split(',') if i != 'NA' and 
                      (i.split(':')[1] == 'NA' or int(i.split(':')[1]) >= MinBaseQuality) and 
                      (line[0],i.split(":")[0]) not in sites_to_ignore and
                      (line[0],i.split(":")[0]) not in ignore_indels]
        hp2s= [(line[0],i.split(":")[0],i.split(':')[2]) for i in line[8].split(',') if i != 'NA' and 
                      (i.split(':')[1] == 'NA' or int(i.split(':')[1]) >= MinBaseQuality) and 
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
        if int(line[6].split(':')[1]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        key= (line[0],line[3],line[4])
        variants= [(line[0],i.split(":")[0],i.split(':')[2]) for i in 
                    set(line[7].split(',') + line[8].split(',') + line[9].split(','))
                    if i != 'NA' and 
                    (i.split(':')[1] == 'NA' or int(i.split(':')[1]) >= MinBaseQuality)]
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
    all_het_snvs= 0
    phased_het_snvs= 0
    pofo_het_snvs= 0
    all_het_indels= 0
    phased_het_indels= 0
    pofo_het_indels= 0
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
                    if len(line[3]) == 1 and len(line[4]) == 1:
                        all_het_snvs += 1
                    else:
                        all_het_indels += 1
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
                        if len(line[3]) == 1 and len(line[4]) == 1:
                            phased_het_snvs += 1
                        else:
                            phased_het_indels += 1
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
                        if len(line[3]) == 1 and len(line[4]) == 1:
                            phased_het_snvs += 1
                        else:
                            phased_het_indels += 1
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
                    if len(line[3]) == 1 and len(line[4]) == 3:
                        all_het_snvs += 1
                    else:
                        all_het_indels += 1
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
                        if len(line[3]) == 1 and len(line[4]) == 3:
                            phased_het_snvs += 1
                        else:
                            phased_het_indels += 1
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
                        if len(line[3]) == 1 and len(line[4]) == 3:
                            phased_het_snvs += 1
                        else:
                            phased_het_indels += 1
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
        if int(line[6].split(':')[1]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        key= (line[0],line[3],line[4])
        variants= [(line[0],str(int(i.split(":")[0])+1),i.split(':')[2]) for i in 
                        set(line[7].split(',') + line[8].split(',') +line[9].split(',')) 
                        if i != 'NA' and 
                        (i.split(':')[1] == 'NA' or int(i.split(':')[1]) >= MinBaseQuality)]
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
    out_freqhp1= out + '_NonPofO_MethylationHP1.tsv'
    out_freqhp2= out + '_NonPofO_MethylationHP2.tsv'
    out_freqhp1_non_pofo = open(out_freqhp1,'w')
    out_freqhp2_non_pofo = open(out_freqhp2,'w')
    freq_header= "Chromosome\tStart\tEnd\tCov\tMod\tFreq"
    out_freqhp1_non_pofo.write(freq_header+"\n")
    out_freqhp2_non_pofo.write(freq_header+"\n")
    p= mp.Pool(threads)
    freq_dicts= p.starmap(out_freq,
                          list(zip(set(list(read_dictHP1.keys())+
                                       list(read_dictHP2.keys())),
                                   repeat(MethylCallfile),
                                   repeat(tool),
                                   repeat(read_dictHP1),
                                   repeat(read_dictHP2),
                                   repeat(callthresh), 
                                   repeat(ref_file))))
    p.close()
    p.join()
    for freq_dict in freq_dicts:
        if freq_dict is None:
            continue
        freq_hp1,freq_hp2 = freq_dict
        if freq_hp1 is not None and freq_hp2 is not None:
            for key,val in freq_hp1.items():
                all_call= len(val)
                mod_call= sum(val)
                out_freqhp1_non_pofo.write('\t'.join(key)+'\t'+str(all_call)+'\t'+
                                            str(mod_call)+'\t'+
                                            str(mod_call/all_call)+'\n')
            freq_hp1.clear()
            for key,val in freq_hp2.items():
                all_call= len(val)
                mod_call= sum(val)
                out_freqhp2_non_pofo.write('\t'.join(key)+'\t'+str(all_call)+'\t'+
                                       str(mod_call)+'\t'+
                                           str(mod_call/all_call)+'\n')
                    
            freq_hp2.clear()
    out_freqhp1_non_pofo.close()
    out_freqhp2_non_pofo.close()
    try:
        subprocess.run("{} {} {} {} {} {} {} {} {} {} {}".format("Rscript",
                                                    os.path.join(os.path.dirname(
                                                                os.path.realpath(__file__)
                                                                    ),
                                                              "DMA_UsingDSS.R"),
                                                    out_freqhp1,
                                                    out_freqhp2,
                                                    out,
                                                    args.equal_disp,
                                                    args.smoothing_flag,
                                                    args.smoothing_span,
                                                    args.delta_cutoff,
                                                    args.pvalue,
                                                    threads),
                       shell=True,
                       check=True)
        subprocess.run("bgzip -f {0} && tabix -f -S 1 -p bed {0}.gz"
                       "".format(out+"_DMLtest.tsv"),
           shell=True,
           check=True)
        subprocess.run("bgzip -f {0} && tabix -f -S 1 -p bed {0}.gz"
                       "".format(out+"_callDML.tsv"),
           shell=True,
           check=True)
    except:
        raise Exception("python subprocess failed."
                        " Make sure you have R, DSS R package, bgzip and"
                        " tabix installed. Moreover, it might be caused because you"
                        " specified both --smoothing_flag and --equal_disp options as FALSE.")
        
    out_pofo = out + '_PofO_Assignment.vcf'
    out_pofo_reads = out + '_PofO_Assignment_reads.tsv'
    assignment_file= open(out_pofo, 'w')
    reads_PofO= open(out_pofo_reads,'w')
    chrom_list= sorted(chrom_list)
    chrom_hp_origin_count= PofO_dmr(known_dmr, 
                                        chrom_list,
                                        out,
                                        min_cg,
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
                    if line[9].startswith('0|1') or line[9].startswith('1|0'):
                        if len(line[3]) == 1 and len(line[4]) == 1:
                            pofo_het_snvs += 1
                        else:
                            pofo_het_indels += 1
                    elif line[9].startswith('1|2'):
                        if len(line[3]) == 1 and len(line[4]) == 3:
                            pofo_het_snvs += 1
                        else:
                            pofo_het_indels += 1
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
    out_scores= open(out + '_PofO_Scores.tsv','w')
    out_scores.write("Chromosome\tOrigin_HP1\tOrigin_HP2\tPofO_Assignment_Score\n")
    out_freqMaternal= out + '_MethylationMaternal.tsv'
    out_freqPaternal= out + '_MethylationPaternal.tsv'
    out_freqMaternal_non_pofo = open(out_freqMaternal,'w')
    out_freqPaternal_non_pofo = open(out_freqPaternal,'w')
    out_freqMaternal_non_pofo.write(freq_header+"\n")
    out_freqPaternal_non_pofo.write(freq_header+"\n")
    for chrom,val in chrom_hp_origin.items():
        for hp,score in val.items():
            if hp=="HP1":
                origin_hp1= score[0]
                if origin_hp1 == "maternal":
                    out_pofo_freq(out_freqhp1,
                                  chrom,
                                  out_freqMaternal_non_pofo)
                    out_pofo_freq(out_freqhp2,
                                  chrom,
                                  out_freqPaternal_non_pofo)
                    origin_hp1= "Maternal"
                    origin_hp2= "Paternal"
                elif origin_hp1 == "paternal":
                    out_pofo_freq(out_freqhp1,
                                  chrom,
                                  out_freqPaternal_non_pofo)
                    out_pofo_freq(out_freqhp2,
                                  chrom,
                                  out_freqMaternal_non_pofo)
                    origin_hp1= "Paternal"
                    origin_hp2= "Maternal"
                out_scores.write('\t'.join([chrom,origin_hp1,origin_hp2,str(score[1]),'\n']))
    out_scores.close()
    out_freqMaternal_non_pofo.close()
    out_freqPaternal_non_pofo.close()
    warnings.warn("Out of {} and {} hetrozygous SNVs and indels in input vcf file"
                  ", {} ({}%) and {} ({}%) could be rephased and "
                  "{} ({}%) and {} ({}%) could be PofO assigned."
                  "".format(all_het_snvs, all_het_indels,
                            phased_het_snvs, round((phased_het_snvs/all_het_snvs)*100,2), 
                            phased_het_indels, round((phased_het_indels/all_het_indels)*100,2),
                            pofo_het_snvs, round((pofo_het_snvs/all_het_snvs)*100,2),
                            pofo_het_indels, round((pofo_het_indels/all_het_indels)*100,2)))
 
    

"""
Specific argument parser.
"""
parser = argparse.ArgumentParser(prog='PatMat.py', add_help=False,
                                 description="Phasing reads and Methylation "
                                             "using strand-seq and nanopore to determine "
                                             "PofO of each homologous chromosome "
                                             "in a single sample.")
required = parser.add_argument_group("Required arguments")
required.add_argument("--bam", "-b",
                      action="store",
                      type=str,
                      required=True,
                      help="The path to the cordinate sorted bam file.")
required.add_argument("--output", "-o",
                      action="store",
                      type=str,
                      required=True,
                      help=("The path to directory and prefix to save "
                            "files. e.g path/to/directory/prefix"))
required.add_argument("--vcf", "-v",
                  action="store",
                  type=str,
                  required=True,
                  default= None,
                  help="The path to the vcf file.")
required.add_argument("--strand_vcf", "-sv",
                      action="store",
                      type=str,
                      required=True,
                      help="The path to the chromosome-scale phased vcf file."
                           "This is the input vcf file that has been phased using strand-seq data.")
required.add_argument("--reference", "-ref",
                      action="store",
                      type=str,
                      required=True,
                      help=("The path to the reference file. File must be indexed"
                            " using samtools faidx."))
required.add_argument("--methylcallfile", "-mc",
                      action="store",
                      type=str,
                      required=True,
                      help=("The path to the per-read methylation "
                            "call file or the bam file with mthylation tag."))
optional = parser.add_argument_group("Optional arguments")
optional.add_argument("--tool_and_callthresh", "-tc",
                           action="store",
                           type=str,
                           required=False,
                           default="guppy:0.4",
                           help=("Software you have used for methylation calling "
                                 "(nanoplish, megalodon, deepsignal,guppy (bam file with 5mC tag)):"
                                 "methylation call threshold for considering a site as "
                                 "methylated, unmethylated or ambiguous in methylation call file. "
                                 "For example, nanopolish:1.5 is when methylation"
                                 " calling performed by nanopolish and a CpG with llr >= 1.5 will be considered "
                                 "as methylated and llr <= -1.5 as unmethylated, anything "
                                 "in between will be considered as ambiguous call and ignored. "
                                 "For megalodon, deepsignl and guppy, call thresold will be delta probability (0-1)"
                                 ". For example threshold 0.6 means >=0.8 is methylated and <=0.2 is not and between"
                                 " 0.2-0.8 will be ignored. Default is guppy:0.4"))
optional.add_argument("--known_dmr", "-kd",
                  action="store",
                  type=str,
                  required= False,
                  default= os.path.join(os.path.dirname(
                                                os.path.realpath(__file__)
                                                    ),
                                             "Imprinted_DMR_List_V1.tsv"),
                  help="The path to the input file for known imprinted DMRs. "
                        "File must have the following information in the following column order: "
                        "chromosome\tstart\tend\tMethylatedAlleleOrigin "
                        "where the methylated allele origin must be either "
                        "maternal or paternal (First row must be header). "
                        "By default, we use iDMR list in repo's patmat directory.")
optional.add_argument("--whatshap_vcf", "-wv",
                      action="store",
                      type=str,
                      required=False,
                      default=None,
                      help=("Path to the WhatsHap phased vcf file that is produced from "
                            "phasing input vcf file using nanopore reads via WhatsHap. This can be useful "
                            "when the chromosome-scale phased variants are very sparce. "
                            "File must be sorted and indexed using tabix."))
optional.add_argument("--whatshap_block", "-wb",
                      action="store",
                      type=str,
                      required=False,
                      default=None,
                      help=("Path to the WhatsHap block file. This file can be "
                            "created using whatshap stats command. File must be " 
                            "converted to a bed format with chromosome\tstart\tend in "
                            "the first three columns (First row must be header). If no block file is given"
                            " then the assumption is that the last part after : sign "
                            "in the 10th column is the phase set (PS) name and blocks will be"
                            " calculated internally."))
optional.add_argument("--black_list", "-bl",
                      action="store",
                      type=str,
                      required=False,
                      default= None,
                      help="List of regions to ignore phased varinats at them."
                      " Three first columns must be chromosome\tstart\tend."
                      " If black list is given the vcf file must be indexed using tabix.")
optional.add_argument("--hapratio", "-hr",
                      action="store",
                      type=float,
                      required=False,
                      default=0.75,
                      help=("0-1 . Minimmum ratio of variants a read must have from a haplotype"
                            " to assign it to that haplotype. Default is 0.75. Note that if you also provide "
                            "WhatsHap phased vcf file this option will be also used to correct phased-block switches"
                            " using Strand-seq phased variants. In this case, it is minimum ratio of phased variants "
                            "at a block that supports the dicision based on strand-seq phased varinats."))
optional.add_argument("--min_base_quality", "-mbq",
                      action="store",
                      type=int,
                      required=False,
                      default=7,
                      help=("Only include bases with phred score higher or"
                            " equal than this option. Default is >=7. if your bam "
                            "does not incude base quality data or cannot be read "
                            "this option will not be use."))
optional.add_argument("--mapping_quality", "-mq",
                      action="store",
                      type=int,
                      required=False,
                      default=20,
                      help=("An integer value to specify thereshold for "
                            "filtering reads based om mapping quality. "
                            "Default is >=20"))
optional.add_argument("--min_variant", "-mv",
                      action="store",
                      type=int,
                      required=False,
                      default=1,
                      help=("minimum number of phased variants must a read "
                            "have to be phased. Default= 1. Note that if you also provide "
                            "WhatsHap phased vcf file this option will be also used to correct phased-block switches"
                            " using Strand-seq phased variants. In this case, it is the minimum number of phased "
                            "variants at a block that need to support the dicision based on strand-seq phased varinats."))
optional.add_argument("--min_read_number", "-mr",
                      action="store",
                      type=int,
                      required=False,
                      default=2,
                      help=("minimum number of reads to support a variant"
                            " to assign to each haplotype. Default= 2"))
optional.add_argument("--min_cg", "-mcg",
                      action="store",
                      type=int,
                      required=False,
                      default= 5,
                      help=("Minimmum number of CpGs an iDMR must have to "
                            " consider it for PofO assignment. Default is 5."))
optional.add_argument("--cpg_difference", "-cd",
                      action="store",
                      type=float,
                      required=False,
                      default= 0.1,
                      help=("Minimum cut off for the fraction of CpGs between haplotypes "
                            "must be differentially methylated at an iDMR to "
                            "consider it for PofO assignment. Default is 0.1."))
optional.add_argument("--threads", "-t",
                      action="store",
                      type=int,
                      required=False,
                      default=4,
                      help="Number of parallel processes. Default is 4.")
optional.add_argument("--chunk_size", "-cs",
                      action="store",
                      type=int,
                      required=False,
                      default=100,
                      help=("Chunk per process. Default is 100"))
optional.add_argument("--include_supplementary", "-is",
                      action="store_true",
                      required=False,
                      help="Also include supplementary reads (Not recommended).")
optional.add_argument("--include_indels", "-ind",
                      action="store_true",
                      required=False,
                      help="Also include indels for read phasing to haplotypes.")

optional.add_argument("--per_read", "-pr",
                      action="store",
                      type=str,
                      required=False,
                      default= None,
                      help="If it is your second try and you have per "
                      "read info file give the path to the per "
                      "read info file. This will be significantly faster."
                      " This is also useful when you want to try different thresholds for options (Note that if you "
                      "also provided WhatsHap phased vcf in your first try, then you cannot use per-read to try "
                      "different --min_variant or --hapratio because these options will be also used to correct"
                      " WhatsHap phased-block switches using strand-seq phased variants),"
                      " different dmr list, black list, include/exclude indels, and include/exclude supp reads.")
optional = parser.add_argument_group("Optional arguments: The following options are DSS options for differential methylation"
                                     " analysis to find differentially methylated CpGs between haplotypes")
optional.add_argument("--delta_cutoff", "-dc",
                            action="store",
                            type=float,
                            default=0.1,
                            required=False,
                            help=("0-1. A threshold for defining differentially "
                                  "methylated loci (DML) or CpGs."
                                  "In DML testing procedure, hypothesis test that the two groups "
                                  "means are equal is conducted at each CpG site. Here if delta is "
                                  "specified, the function will compute the posterior probability that "
                                  "the difference of the means are greater than delta,"
                                  " and then call DML based on that. Default is 0.1."))
optional.add_argument("--pvalue", "-pv",
                      action="store",
                      type=float,
                      required=False,
                      default= 0.001,
                      help=("0-1. When delta is not specified, this is the threshold of p-value for defining DML and "
                            "loci with p-value less than this threshold will be deemed DMLs."
                            " When delta is specified, CpG sites with posterior probability greater than 1-pvalue"
                            "_threshold are deemed DML. Default is 0.001"))
optional.add_argument("--smoothing_span", "-sms",
                            action="store",
                            type=int,
                            default=500,
                            required=False,
                            help=("The size of smoothing window, in "
                                  "basepairs. Default is 500."))
optional.add_argument("--smoothing_flag", "-smf",
                        action="store",
                        type=str,
                        default="TRUE",
                        required=False,
                        help=("TRUE/FALSE. A flag to indicate whether to apply smoothing"
                              " in estimating mean methylation levels."
                              " For more instruction see DSS R package guide. Default is TRUE."))
optional.add_argument("--equal_disp", "-ed",
                        action="store",
                        type=str,
                        default="FALSE",
                        required=False,
                        help=("TRUE/FALSE. A flag to indicate whether the "
                              "dispersion in two groups are deemed equal or not. "
                              "For more instruction see DSS R package guide. Default is FALSE."
                              " Because there is no biological replicate here, you should"
                              " specify either equal_disp TRUE or smoothing_flag TRUE. Do not specify both as FALSE."))
optional = parser.add_argument_group("Help and version options")
optional.add_argument('--version', action='version', 
                      version='%(prog)s 1.3.0_dev',
                      help= "Print program's version and exit")
optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                      help="Print this help and exit.")
args = parser.parse_args()

if __name__ == '__main__':
    main(args)
