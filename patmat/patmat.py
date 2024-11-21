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
PatMat: Parent-of-origin (Paternal and Maternal) resolved chromosome-scale haplotyping
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
from tqdm import tqdm

def get_variant_info(feed_list,
                     alignment_file,
                     chrom):
    '''
    This function maps each read to heterozygous variants and returns a list of  
    haplotype 1, haplotype 2, and unphased variants mapped to each read.
    '''
    read_var_list= list()
    read_hap_list= list()
    read_var_cov= defaultdict(list)
    samfile = pysam.AlignmentFile(alignment_file, 'rb')
    for varinfo in feed_list:
        gt,position,ref,alt= varinfo
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
                    read_mq = pileupread.alignment.mapping_quality
                    if read_mq < args.mapping_quality:
                        continue
                    if pileupread.alignment.is_supplementary and not args.include_supplementary:
                        continue
                    read_id= pileupread.alignment.query_name
                    flag= pileupread.alignment.flag
                    read_start= pileupread.alignment.reference_start
                    ext_key= (read_id,str(flag),str(read_mq),str(read_start+1))
                    read_var_cov[(chrom,position)].append((chrom,read_id))
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    read_base= None
                    read_base1= None
                    read_base2= None
                    if gt != "1/2":                    
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
                        if read_base == ref:
                            read_var_list.append((chrom,position,read_base.upper(),*ext_key))
                            if gt == '1|0':
                                read_hap_list.append((chrom,*ext_key,2))
                            elif gt == '0|1' and read_base == ref:
                                read_hap_list.append((chrom,*ext_key,1))
                        elif read_base == alt:
                            read_var_list.append((chrom,position,read_base.upper(),*ext_key))
                            if gt == '1|0' and read_base == alt:
                                read_hap_list.append((chrom,*ext_key,1))
                            elif gt == '0|1' and read_base == alt:
                                read_hap_list.append((chrom,*ext_key,2))
                    else:
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
                        if read_base1 == alt1 or read_base1 == alt2:
                            read_var_list.append((chrom,position,read_base1.upper(),*ext_key))
                        elif read_base2 == alt1 or read_base2 == alt2:
                            read_var_list.append((chrom,position,read_base2.upper(),*ext_key))
                            
    return read_var_list, read_hap_list, read_var_cov



def getChromsFromBAM(filename):
    chroms = set()
    stats = pysam.idxstats(filename)
    for row in stats.split("\n"):
        fields = row.split("\t")
        if fields[0] != '*' and fields[0] != '':
            chroms.add(fields[0])
    return chroms



def get_sv_pofo(feed_list,
                   alignment_file):
    '''
    This function maps each read to heterozygous variants and returns a list of  
    haplotype 1, haplotype 2, and unphased variants mapped to each read.
    '''
    sv_hp_dict1= defaultdict(int)
    sv_hp_dict2= defaultdict(int)
    samfile = pysam.AlignmentFile(alignment_file, 'rb')
    for varinfo in feed_list:
        chrom,position,sv_reads = varinfo
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
                    read_id= pileupread.alignment.query_name
                    try:
                        hp_tag= pileupread.alignment.get_tag("HP")
                    except:
                        continue
                    if read_id in sv_reads:
                        if hp_tag==1:
                            sv_hp_dict1[(chrom,position)] += 1
                            sv_hp_dict2[(chrom,position)] += 0
                        elif hp_tag==2:
                            sv_hp_dict1[(chrom,position)] += 0
                            sv_hp_dict2[(chrom,position)] += 1
                
    return sv_hp_dict1, sv_hp_dict2



def openalignment(alignment_file,
                  window):
    '''
    Opens an alignment file and creates bam iterator
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
             out,
             out_meth,
             min_cg,
             cpg_difference):
    '''
    This function maps differentially methylated CpGs 
    to the known iDMRs for PofO assignment 
    '''
    tb_calldml= tabix.open(out+"_Temp_callDML.tsv.gz")
    tb_dmltest= tabix.open(out+"_Temp_DMLtest.tsv.gz")
    dmr_file= openfile(known_dmr)
    header= next(dmr_file).rstrip()
    chrom_hp_origin_count= defaultdict(lambda: defaultdict(int))
    chrom_hp_origin_count_all_cg= defaultdict(lambda: defaultdict(int))
    chrom_hp_origin_count_diff_cg= defaultdict(lambda: defaultdict(int))
    chrom_hp_origin_count_dmr= defaultdict(lambda: defaultdict(int))
    for line in dmr_file:
        line= line.rstrip().split('\t')
        dmr_chrom= line[0]
        dmr_start= int(line[1]) -1
        dmr_end= int(line[2])
        origin= line[3].lower()
        if origin not in ["maternal","paternal"]:
            warnings.warn("iDMR: {} does not have a valid origin "
                          "(must be paternal or maternal. case insensitive). "
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
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+
                           '\tNA\tNA\tNA\tNA\tIgnored:No_CpG\n')
            continue
        hp1_freq= hp1_freq/num_cg
        hp2_freq= hp2_freq/num_cg
        if num_cg < min_cg and abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+
                           str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + 
                           'Ignored:DidNotMeet_min_cg_And_cpg_difference\n')
            continue
        elif num_cg < min_cg and diff_cg_hp1 == 0 and diff_cg_hp2 == 0:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+
                           str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + 
                           'Ignored:DidNotMeet_min_cg_And_DifferentiallyMethylated'
                           'CpGsIsZeroInBothHaplotypes\n')
            continue
        elif num_cg < min_cg:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+
                           str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' +
                           'Ignored:DidNotMeet_min_cg\n')
            continue
        elif abs(diff_cg_hp1 - diff_cg_hp2) / num_cg < cpg_difference:
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+
                           str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' +
                           'Ignored:DidNotMeet_cpg_difference\n')
            continue
        elif diff_cg_hp1 == 0 and diff_cg_hp2 == 0: 
            out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+
                           str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+ 
                           '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + 
                           'Ignored:DifferentiallyMethylatedCpGsIsZeroInBothHaplotypes\n')
            continue
        out_meth.write('\t'.join(line)+'\t'+str(num_cg)+'\t'+
                       str(diff_cg_hp1)+'\t'+str(diff_cg_hp2)+
                       '\t'+str(hp1_freq) + '\t' + str(hp2_freq)+'\t' + 
                       'Included\n')  
        if origin == 'maternal':
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += diff_cg_hp1
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += diff_cg_hp2
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += diff_cg_hp1
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += diff_cg_hp2
            if diff_cg_hp1 > diff_cg_hp2:
                chrom_hp_origin_count_dmr[(dmr_chrom, 'maternal')]['HP1'] += 1
                chrom_hp_origin_count_dmr[(dmr_chrom, 'paternal')]['HP2'] += 1
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'maternal')]['HP1'] += num_cg
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'paternal')]['HP2'] += num_cg
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'maternal')]['HP1'] += diff_cg_hp1+diff_cg_hp2
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'paternal')]['HP2'] += diff_cg_hp1+diff_cg_hp2
            elif diff_cg_hp1 < diff_cg_hp2:
                chrom_hp_origin_count_dmr[(dmr_chrom, 'maternal')]['HP2'] += 1
                chrom_hp_origin_count_dmr[(dmr_chrom, 'paternal')]['HP1'] += 1
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'maternal')]['HP2'] += num_cg
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'paternal')]['HP1'] += num_cg
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'maternal')]['HP2'] += diff_cg_hp1+diff_cg_hp2
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'paternal')]['HP1'] += diff_cg_hp1+diff_cg_hp2
        if origin == 'paternal':
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP1'] += diff_cg_hp2 
            chrom_hp_origin_count[(dmr_chrom, 'maternal')]['HP2'] += diff_cg_hp1
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP2'] += diff_cg_hp2 
            chrom_hp_origin_count[(dmr_chrom, 'paternal')]['HP1'] += diff_cg_hp1
            if diff_cg_hp1 > diff_cg_hp2:
                chrom_hp_origin_count_dmr[(dmr_chrom, 'maternal')]['HP2'] += 1
                chrom_hp_origin_count_dmr[(dmr_chrom, 'paternal')]['HP1'] += 1
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'maternal')]['HP2'] += num_cg
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'paternal')]['HP1'] += num_cg
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'maternal')]['HP2'] += diff_cg_hp1+diff_cg_hp2
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'paternal')]['HP1'] += diff_cg_hp1+diff_cg_hp2
            elif diff_cg_hp1 < diff_cg_hp2:
                chrom_hp_origin_count_dmr[(dmr_chrom, 'maternal')]['HP1'] += 1
                chrom_hp_origin_count_dmr[(dmr_chrom, 'paternal')]['HP2'] += 1
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'maternal')]['HP1'] += num_cg
                chrom_hp_origin_count_all_cg[(dmr_chrom, 'paternal')]['HP2'] += num_cg
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'maternal')]['HP1'] += diff_cg_hp1+diff_cg_hp2
                chrom_hp_origin_count_diff_cg[(dmr_chrom, 'paternal')]['HP2'] += diff_cg_hp1+diff_cg_hp2
    dmr_file.close()
    out_meth.close()
    return (chrom_hp_origin_count, 
            chrom_hp_origin_count_all_cg, 
            chrom_hp_origin_count_diff_cg,
            chrom_hp_origin_count_dmr)


def out_freq_methbam(out,
                     processes,
                     reference,
                     pbcg,
                     pb_tech):
    out_freqhp1= out + '_Temp_NonPofO_HP1-HP2_MethylationHP1.tsv'
    out_freqhp2= out + '_Temp_NonPofO_HP1-HP2_MethylationHP2.tsv'
    out_dir=os.path.dirname(out)
    out_pref= os.path.basename(out)
    if not pb_tech:
        subprocess.run("modkit pileup -t {} --prefix {} "
                       "--partition-tag HP --combine-strands --cpg -r {} "
                       "{} {}".format(processes,out_pref+"_Temp-NonPofO_CpGModFreq",
                                      reference, out+"_Temp-NonPofO_dmr.bam",
                                      out_dir),
                       shell=True,
                       check=True)
        
        subprocess.run("awk -F'\t' '$4==\"m\" {{print $1,$2,$3,$10,$12,$12/$10}}' OFS='\t' {} | "
                       "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
                        "".format(out+"_Temp-NonPofO_CpGModFreq_1.bed",
                                          out_freqhp1),
                        shell=True,
                        check=True)
        subprocess.run("awk -F'\t' '$4==\"m\" {{print $1,$2,$3,$10,$12,$12/$10}}' OFS='\t' {} | "
                       "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
                        "".format(out+"_Temp-NonPofO_CpGModFreq_2.bed",
                                          out_freqhp2),
                        shell=True,
                        check=True)
    else:
        subprocess.run("aligned_bam_to_cpg_scores --bam {} --output-prefix {}"
                       " --model {} --threads {} --modsites-mode reference "
                       "--ref {}".format(out+"_Temp-NonPofO_dmr.bam",
                                         out+"_Temp-NonPofO_CpGModFreq",
                                         pbcg,
                                         processes,reference),
                       shell=True,
                       check=True)
        
        subprocess.run("awk -F'\t' '{{print $1,$2,$3,$6,$7,$7/$6}}' OFS='\t' {} | "
                       "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
                        "".format(out+"_Temp-NonPofO_CpGModFreq.hap1.bed",
                                          out_freqhp1),
                        shell=True,
                        check=True)
        subprocess.run("awk -F'\t' '{{print $1,$2,$3,$6,$7,$7/$6}}' OFS='\t' {} | "
                       "sed '1i Chromosome\tStart\tEnd\tCov\tMod\tFreq' > {} "
                        "".format(out+"_Temp-NonPofO_CpGModFreq.hap2.bed",
                                          out_freqhp2),
                        shell=True,
                        check=True)
    subprocess.run("rm {}*".format(out+"_Temp-NonPofO_CpGModFreq"),
                    shell=True,
                    check=True)
                  

def out_pofo_freq(hp_fre,
                  chrom,
                  output):
    '''
    Outputs PofO assigned methylation frequency files.
    '''
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


def vcf2dict(vcf_strand,
             chrom):
    """
    Process and converts the strand-seq vcf file to a dictionary
    for downstream use.
    """
    strand_phased_vars= 0
    vcf_dict= defaultdict(dict)
    vcf_file = openfile(vcf_strand)
    for line in vcf_file:
        line = line.rstrip().split('\t')
        if line[0]!=chrom:
            continue
        if line[9].startswith(('.|0','0|.','1|.','.|1')):
            line[9] = line[9].replace(".|0","1|0"
                                                ).replace(".|1","0|1"
                                                          ).replace("1|.","1|0"
                                                                    ).replace("0|.","0|1")
            warnings.warn("{}:{} variant in strand-seq vcf has .|0 or 0|. "
                          "or .|1 or 1|. phased genotype. Note that it "
                          "will be interpreted as 1|0 or 0|1 or 0|1 or "
                          "1|0".format(line[0], line[1]))
        if line[9].startswith(('1|0','0|1')):
            strand_phased_vars += 1
            vcf_dict[line[0]][line[1]] = line[9].split(':')[0]
    vcf_file.close()
    return vcf_dict, strand_phased_vars


def alignment_writer(bam,
                     chrom,
                     reads_hap,
                     chrom_hp_origin,
                     outfile):
    bamiter = bam.fetch(chrom)
    if chrom in chrom_hp_origin and chrom_hp_origin[chrom]['HP1'][0] == 'maternal':
        for read in bamiter:
            if (read.mapping_quality < args.mapping_quality or 
                read.is_secondary or read.is_qcfail or 
                read.is_duplicate or read.is_unmapped or
                (read.is_supplementary and not args.include_supplementary)):
                outfile.write(read)
                continue
            read.set_tag("HP", None)
            read.set_tag("PS", None)
            read_id = read.query_name
            ref_name = read.reference_name
            if(ref_name,read_id) in reads_hap:
                read.set_tag("HP", int(reads_hap[(ref_name,read_id)]))
            outfile.write(read)
    elif chrom in chrom_hp_origin and chrom_hp_origin[chrom]['HP2'][0] == 'maternal':
        for read in bamiter:
            if (read.mapping_quality < args.mapping_quality or 
                read.is_secondary or read.is_qcfail or 
                read.is_duplicate or read.is_unmapped or
                (read.is_supplementary and not args.include_supplementary)):
                outfile.write(read)
                continue
            read.set_tag("HP", None)
            read.set_tag("PS", None)
            read_id = read.query_name
            ref_name = read.reference_name
            if(ref_name,read_id) in reads_hap:
                if reads_hap[(ref_name,read_id)]==1:
                    read.set_tag("HP", 2)
                else:
                    read.set_tag("HP", 1)
            outfile.write(read)
    else:
        for read in bamiter:
            if (read.mapping_quality < args.mapping_quality or 
                read.is_secondary or read.is_qcfail or 
                read.is_duplicate or read.is_unmapped or
                (read.is_supplementary and not args.include_supplementary)):
                outfile.write(read)
                continue
            read.set_tag("HP", None)
            read.set_tag("PS", None)
            outfile.write(read)

   
    
def pofo_sv_write(sv_file,
                  out,
                  chrom_hp_origin,
                  reads_hap,
                  min_read_reassignment):
    sv_assignment_file= open(out +"_"+os.path.basename(sv_file) + 
                                      '_PofO_Assignment_SVs.vcf', 'w')
    sv_assignment_file_info= open(out +"_"+os.path.basename(sv_file) + 
                                      '_PofO_Assignment_SVs_info.tsv', 'w') 
    with openfile(sv_file) as vf:
        for line in vf:
            if line.startswith("##"):
                sv_assignment_file.write(line)
                continue
            elif line.startswith("#"):
                sv_assignment_file.write(line)
                sv_assignment_file_info.write(line.rstrip()+"\tNumHp1ReadsFromColumn8\t"
                                         "NumHp2ReadsFromColumn8"
                                         "\tNumMaternalReadsFromColumn8"
                                         "\tNumPaternalReadsFromColumn8\n")
                continue
            line=line.rstrip().split("\t")
            gt=line[9].split(":")[0]
            if 'PS' in line[8].split(":"):
                ps_index= line[8].split(":").index('PS')
                new_ps= line[8].split(":")
                new_hp= line[9].split(":")
                new_ps.pop(ps_index)
                new_hp.pop(ps_index)
            else:
                new_ps= line[8].split(":")
                new_hp= line[9].split(":")
            if not args.include_all_variants and line[6] not in ["PASS","."]:
                line_out= (line[0:8]+
                            [':'.join(new_ps)]+
                            [':'.join(new_hp).replace("|", "/").
                                              replace("1/0", "0/1").
                                              replace("2/1", "1/2")]+
                            ["NA"]*4)
                sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                sv_assignment_file_info.write('\t'.join(line_out)+'\n')
                continue
            if line[0] in chrom_hp_origin:
                if gt in ["0/1","1/0","0|1","1|0"] and "RNAMES=" in line[7]:
                    hp1_count= 0
                    hp2_count= 0
                    for read_ID in line[7].split("RNAMES=")[1].split(";")[0].split(","):
                        if (line[0],read_ID) in reads_hap:
                            if reads_hap[(line[0],read_ID)]==1:
                                hp1_count+= 1
                            elif reads_hap[(line[0],read_ID)]==2:
                                hp2_count+= 1
                    if (hp1_count > hp2_count and
                        hp1_count >= min_read_reassignment):
                        if chrom_hp_origin[line[0]]['HP1'][0] == 'maternal':
                            line_out= (line[0:8]+
                                        [':'.join(new_ps)+":PS"]+
                                        ["1|0:"+':'.join(new_hp[1:])+":Pat|Mat"]+
                                        [str(hp1_count), str(hp2_count),
                                          str(hp2_count), str(hp1_count)])
                            sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                            sv_assignment_file_info.write('\t'.join(line_out)+'\n')
                        if chrom_hp_origin[line[0]]['HP1'][0] == 'paternal':
                            line_out= (line[0:8]+
                                        [':'.join(new_ps)+":PS"]+
                                        ["0|1:"+':'.join(new_hp[1:])+":Mat|Pat"]+
                                        [str(hp1_count), str(hp2_count),
                                          str(hp1_count), str(hp2_count)])
                            sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                            sv_assignment_file_info.write('\t'.join(line_out)+'\n')
                    elif (hp2_count > hp1_count and
                          hp2_count >= min_read_reassignment):
                        if chrom_hp_origin[line[0]]['HP2'][0] == 'maternal':
                            line_out= (line[0:8]+
                                        [':'.join(new_ps)+":PS"]+
                                        ["1|0:"+':'.join(new_hp[1:])+":Pat|Mat"]+
                                        [str(hp1_count), str(hp2_count),
                                          str(hp1_count), str(hp2_count)])
                            sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                            sv_assignment_file_info.write('\t'.join(line_out)+'\n')
                        if chrom_hp_origin[line[0]]['HP2'][0] == 'paternal':
                            line_out= (line[0:8]+
                                        [':'.join(new_ps)+":PS"]+
                                        ["0|1:"+':'.join(new_hp[1:])+":Mat|Pat"]+
                                        [str(hp1_count), str(hp2_count),
                                          str(hp2_count), str(hp1_count)])
                            sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                            sv_assignment_file_info.write('\t'.join(line_out)+'\n')
                    else:
                        line_out= (line[0:8]+
                                    [':'.join(new_ps)]+
                                    [line[9].replace("|", "/").
                                              replace("1/0", "0/1")]+
                                    [str(hp1_count), str(hp2_count),
                                      "NA","NA"])
                        sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                        sv_assignment_file_info.write('\t'.join(line_out)+'\n')                                    
                else:
                    line_out= (line[0:8]+
                                [':'.join(new_ps)]+
                                [line[9].replace("|", "/").
                                  replace("1/0", "0/1").
                                  replace("2/1", "1/2")]+
                                ["NA"]*4)
                    sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                    sv_assignment_file_info.write('\t'.join(line_out)+'\n')
            else:
                line_out= (line[0:8]+
                            [':'.join(new_ps)]+
                            [line[9].replace("|", "/").
                              replace("1/0", "0/1").
                              replace("2/1", "1/2")]+
                            ["NA"]*4)
                sv_assignment_file.write('\t'.join(line_out[0:-4])+'\n')
                sv_assignment_file_info.write('\t'.join(line_out)+'\n')
    sv_assignment_file.close()
    sv_assignment_file_info.close()
    
    
def strand_vcf2dict_phased(vcf_strand,
                           vcf,
                           include_all_variants,
                           chrom):
    '''
    Intersects input vcf and strand-seq vcf and stores phased 
    and unphased heterozygous variants into a dictionary for read phasing.
    '''
    final_dict= defaultdict(set)
    vcf_dict, strand_phased_vars= vcf2dict(vcf_strand,chrom)
    with openfile(vcf) as vf:
        for line in vf:
            line=line.rstrip().split('\t')
            if line[0]!=chrom:
                continue
            if not include_all_variants and line[6] not in ["PASS","."]:
                continue
            if (line[9].startswith(('0/1','1/0','0|1','1|0'))
                and
                (line[0] in vcf_dict and 
                 line[1] in vcf_dict[line[0]])):
                final_dict[line[0]].add((vcf_dict[line[0]][line[1]],
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif line[9].startswith(('0/1','1/0','0|1','1|0')):
                final_dict[line[0]].add(("0/1",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif line[9].startswith(('1/2','1|2','2/1','2|1')):
                final_dict[line[0]].add(("1/2",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        (line[4].split(',')[0].upper(),
                                        line[4].split(',')[1].upper())))
    return final_dict, strand_phased_vars
                

def vcf2dict_phased(blocks, 
                    vcf_strand,
                    hapRatio,
                    minvariantblock,
                    vcf,
                    include_all_variants,
                    chrom):
    '''
    In case --phased option is given, this function uses phased variants
    from strand-seq to correct PofO phasing switches across phased blocks
    . This function then intersects input vcf and strand-seq vcf and
    stores phased and unphased heterozygous variants into a 
    dictionary for read phasing.
    '''
    final_dict= defaultdict(set)
    phase_block_stat= defaultdict(int)
    vcf_dict, strand_phased_vars= vcf2dict(vcf_strand,chrom)
    tb_vcf= tabix.open(os.path.abspath(vcf))
    for block in blocks:
        b_chrom, b_start, b_end= block
        if b_chrom!=chrom:
            continue
        agreement_count= 0
        disagreement_count= 0
        try:
            records_whats = tb_vcf.query(b_chrom, b_start-1, b_end+1)
            records_whats = list(records_whats)
        except:
            warnings.warn("{}:{}-{} block cannot be extracted from vcf file. "
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
            if agreement_count > disagreement_count:
                phase_block_stat[(b_chrom, str(b_start),"agreement")]+=agreement_count
                phase_block_stat[(b_chrom, str(b_start),"disagreement")]+=disagreement_count
            else:
                phase_block_stat[(b_chrom, str(b_start),"agreement")]+=disagreement_count
                phase_block_stat[(b_chrom, str(b_start),"disagreement")]+=agreement_count
            
            for vcf_line in records_whats:
                if vcf_line[1] in vcf_dict[vcf_line[0]]:
                    continue
                if vcf_line[9].startswith(('1|0','0|1')):
                    if (agreement_count > disagreement_count and
                        agreement_count >= minvariantblock and 
                        agreement_count/(agreement_count+disagreement_count) >= hapRatio and
                        vcf_line[1] not in disagreement_sites):
                        vcf_dict[vcf_line[0]][vcf_line[1]] = vcf_line[9].split(':')[0]
                    elif (disagreement_count > agreement_count and
                          disagreement_count >= minvariantblock and 
                          disagreement_count/(agreement_count+disagreement_count) >= hapRatio and
                          vcf_line[1] not in agreement_sites):
                        vcf_dict[vcf_line[0]][vcf_line[1]] = vcf_line[9].split(':')[0][::-1]
    with openfile(vcf) as vf:
        for line in vf:
            line=line.rstrip().split('\t')
            if line[0]!=chrom:
                continue
            if not include_all_variants and line[6] not in ["PASS","."]:
                continue
            if (line[9].startswith(('0/1','1/0','0|1','1|0'))
                and
                (line[0] in vcf_dict and 
                 line[1] in vcf_dict[line[0]])):
                final_dict[line[0]].add((vcf_dict[line[0]][line[1]],
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif line[9].startswith(('0/1','1/0','0|1','1|0')):
                final_dict[line[0]].add(("0/1",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        line[4].upper()))
            elif line[9].startswith(('1/2','1|2','2/1','2|1')):
                final_dict[line[0]].add(("1/2",
                                        int(line[1])-1,
                                        line[3].upper(),
                                        (line[4].split(',')[0].upper(),
                                        line[4].split(',')[1].upper())))

    vcf_dict.clear()
    return final_dict, strand_phased_vars, phase_block_stat


def write_sam_phase(in_file,
                    out_file,
                    reads_hap):
    out_nonpofo_bam= gzip.open(out_file, 'wb')
    with openfile(in_file) as sf:
        for line in sf:
            if line.startswith(("@HD","@SQ","@RG")):
                out_nonpofo_bam.write(line.encode())
                continue
            elif line.startswith(("@")):
                continue
            line= line.rstrip().split('\t')
            out_read= '\t'.join(line)+'\n'
            if (int(line[4]) < args.mapping_quality or 
                line[1] not in ["0","16","2048","2064"]):
                continue
            elif not args.include_supplementary and line[1] in ["2048","2064"]:
                continue
            elif (line[2],line[0]) in reads_hap:
                out_read= '\t'.join(line)+'\t'+"HP:i:"+str(reads_hap[(line[2],line[0])])+'\n'
                out_nonpofo_bam.write(out_read.encode())
    out_nonpofo_bam.close()
    
    

def per_read_variant(vcf_dict,
                     bam_file,
                     chunk,
                     processes):
    '''
    This function extracts per-read information for variants.
    '''
    variant_dict= defaultdict(set)
    read_dict_HP1= defaultdict(int)
    read_dict_HP2= defaultdict(int)
    per_var_cov = defaultdict(list)
    for chrom,feed_list in vcf_dict.items():
        bamiter, bam, count = openalignment(bam_file, chrom)
        if count > 0:
            feed_list = [list(feed_list)[x:x+chunk]
                                  for x in range(0, len(feed_list),
                                                chunk)]
            feed_list = [feed_list[x:x+processes]
                                  for x in range(0, len(feed_list),
                                                processes)]
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
                            read_var_list, read_hap_list, read_var_cov= result
                            for var_info in read_var_list:
                                variant_dict[var_info[0:-4]].add(var_info[-4:])
                            for read_hap in read_hap_list:
                                if read_hap[-1] == 1:
                                    read_dict_HP1[read_hap[0:-1]] += 1
                                    read_dict_HP2[read_hap[0:-1]] += 0
                                elif read_hap[-1] == 2:
                                    read_dict_HP1[read_hap[0:-1]] += 0
                                    read_dict_HP2[read_hap[0:-1]] += 1
                            for key,val in read_var_cov.items():
                                per_var_cov[key] += val
                            
                    pbar.update(1)
        else:
            warnings.warn("{} does not have any mapped reads in alignment "
                          "file Or alignment is truncated or corrupt indexed. "
                          "Skipping it.".format(chrom))
            
    return (variant_dict, read_dict_HP1, 
            read_dict_HP2, per_var_cov)


def pofo_final_dict(chrom_hp_origin_count,
                    chrom_hp_origin_count_dmr,
                    chrom_hp_origin_count_all_cg,
                    chrom_hp_origin_count_diff_cg):
    chrom_hp_origin = defaultdict(dict)
    for key, val in chrom_hp_origin_count.items():
        chrom, origin= key
        if origin.lower() == 'maternal':
            hp1_cg_count= val['HP1']
            hp2_cg_count= val['HP2']
            hp1_dmr_count= chrom_hp_origin_count_dmr[(chrom, 'maternal')]['HP1']
            hp2_dmr_count= chrom_hp_origin_count_dmr[(chrom, 'maternal')]['HP2']
            hp1_allcg_count= chrom_hp_origin_count_all_cg[(chrom, 'maternal')]['HP1']
            hp2_allcg_count= chrom_hp_origin_count_all_cg[(chrom, 'maternal')]['HP2']
            hp1_alldiffcg_count= chrom_hp_origin_count_diff_cg[(chrom, 'maternal')]['HP1']
            hp2_alldiffcg_count= chrom_hp_origin_count_diff_cg[(chrom, 'maternal')]['HP2']
            if hp1_cg_count > hp2_cg_count:
                chrom_hp_origin[chrom]['HP1'] = ['maternal',
                                                  hp1_cg_count, hp2_cg_count,
                                                  hp1_dmr_count, hp2_dmr_count,
                                                  hp1_alldiffcg_count, hp2_alldiffcg_count,
                                                  hp1_allcg_count, hp2_allcg_count]
                                                  
                chrom_hp_origin[chrom]['HP2'] = ['paternal',
                                                  hp1_cg_count, hp2_cg_count,
                                                  hp1_dmr_count, hp2_dmr_count,
                                                  hp1_alldiffcg_count, hp2_alldiffcg_count,
                                                  hp1_allcg_count, hp2_allcg_count]
            elif hp2_cg_count > hp1_cg_count:
                chrom_hp_origin[chrom]['HP2'] = ['maternal',
                                                  hp2_cg_count, hp1_cg_count,
                                                  hp2_dmr_count, hp1_dmr_count,
                                                  hp2_alldiffcg_count, hp1_alldiffcg_count,
                                                  hp2_allcg_count, hp1_allcg_count]
                chrom_hp_origin[chrom]['HP1'] = ['paternal',
                                                  hp2_cg_count, hp1_cg_count,
                                                  hp2_dmr_count, hp1_dmr_count,
                                                  hp2_alldiffcg_count, hp1_alldiffcg_count,
                                                  hp2_allcg_count, hp1_allcg_count]
        elif origin.lower() == 'paternal':
            hp1_cg_count= val['HP1']
            hp2_cg_count= val['HP2']
            hp1_dmr_count= chrom_hp_origin_count_dmr[(chrom, 'paternal')]['HP1']
            hp2_dmr_count= chrom_hp_origin_count_dmr[(chrom, 'paternal')]['HP2']
            hp1_allcg_count= chrom_hp_origin_count_all_cg[(chrom, 'paternal')]['HP1']
            hp2_allcg_count= chrom_hp_origin_count_all_cg[(chrom, 'paternal')]['HP2']
            hp1_alldiffcg_count= chrom_hp_origin_count_diff_cg[(chrom, 'paternal')]['HP1']
            hp2_alldiffcg_count= chrom_hp_origin_count_diff_cg[(chrom, 'paternal')]['HP2']
            if hp1_cg_count > hp2_cg_count:
                chrom_hp_origin[chrom]['HP1'] = ['paternal',
                                                  hp1_cg_count, hp2_cg_count,
                                                  hp1_dmr_count, hp2_dmr_count,
                                                  hp1_alldiffcg_count, hp2_alldiffcg_count,
                                                  hp1_allcg_count, hp2_allcg_count]
                chrom_hp_origin[chrom]['HP2'] = ['maternal',
                                                  hp1_cg_count, hp2_cg_count,
                                                  hp1_dmr_count, hp2_dmr_count,
                                                  hp1_alldiffcg_count, hp2_alldiffcg_count,
                                                  hp1_allcg_count, hp2_allcg_count]
            elif hp2_cg_count > hp1_cg_count:
                chrom_hp_origin[chrom]['HP2'] = ['paternal',
                                                  hp2_cg_count, hp1_cg_count,
                                                  hp2_dmr_count, hp1_dmr_count,
                                                  hp2_alldiffcg_count, hp1_alldiffcg_count,
                                                  hp2_allcg_count, hp1_allcg_count]
                chrom_hp_origin[chrom]['HP1'] = ['maternal',
                                                  hp2_cg_count, hp1_cg_count,
                                                  hp2_dmr_count, hp1_dmr_count,
                                                  hp2_alldiffcg_count, hp1_alldiffcg_count,
                                                  hp2_allcg_count, hp1_allcg_count]
    return chrom_hp_origin





def get_block(vcf):
    '''
    In case --phased option is provided this function extracts "
    "phased blocks from vcf file.
    '''
    blocks_dict= defaultdict(set)
    final= list()
    with openfile(vcf) as vp:
        for line in vp:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            if "|" in line[9].split(":")[0]:
               blocks_dict[(line[0],line[9].split(":")[-1])].add(int(line[1]))
    for key, val in blocks_dict.items():
        val= sorted(val)
        final.append((key[0], val[0], val[-1]))
        blocks_dict[key]= (str(val[0]), str(val[-1]))
    return final, blocks_dict

def main(args):
    '''
    The main function that uses user's inputs and other functions to phase and
    assign PofO to variants and methylation.
    '''
    hapRatio = args.hapratio
    minvariant= args.min_variant
    minvariantblock= args.min_variant_block
    min_cg= args.min_cg
    bam_file = os.path.abspath(args.bam)
    processes = args.processes
    dss_processes= args.dss_processes
    if dss_processes is None:
        dss_processes= processes
    cpg_difference= args.cpg_difference
    chunk = args.chunk_size
    min_read_reassignment= args.min_read_number
    out = os.path.abspath(args.output)
    reference= os.path.abspath(args.reference)
    vcf= os.path.abspath(args.vcf)
    known_dmr= os.path.abspath(args.known_dmr)
    if not os.path.isfile(vcf+".tbi"):
        raise Exception("It seems that vcf file "
                        "is not index by tabix.")
    vcf_tb= tabix.open(vcf)
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
    re_assignment_vars= dict()
    chroms= dict()
    with openfile(vcf) as vf_file:
        for line in vf_file:
            if line.startswith("#"):
                continue
            else:
                line=line.rstrip().split('\t')
                chroms[line[0]]= int(line[1])
    
    bam_choms= getChromsFromBAM(bam_file)
    reads_hap= dict()
    for chrom in sorted(chroms.keys()):
        print("#############  Processing chromosome {}  #############".format(chrom))
        if args.strand_vcf is not None and not args.phased:
            vcf_strand = os.path.abspath(args.strand_vcf)
            final_dict, strand_phased_vars= strand_vcf2dict_phased(vcf_strand, 
                                                           vcf, 
                                                           args.include_all_variants,
                                                           chrom)
        elif args.strand_vcf is not None and args.phased:
            vcf_strand = os.path.abspath(args.strand_vcf)
            blocks_phased, blocks_dict= get_block(vcf)
            (final_dict,
             strand_phased_vars,
             phase_block_stat)= vcf2dict_phased(blocks_phased, 
                                                vcf_strand,
                                                hapRatio,
                                                minvariantblock,
                                                vcf,
                                                args.include_all_variants,
                                                chrom)
        else:
            raise Exception("No strand-seq vcf is given.")
        
        (variant_dict,
         read_dict_HP1, 
         read_dict_HP2,
         per_var_cov)= per_read_variant(final_dict,
                                        bam_file,
                                        chunk,
                                        processes)
        if not read_dict_HP1 and not read_dict_HP1:
            continue
        reads_hap_temp= dict()
        read_hap_reassign= defaultdict(lambda: defaultdict(int))
        variant_dict_HP1_ave= defaultdict(lambda: defaultdict(int))
        variant_dict_HP2_ave= defaultdict(lambda: defaultdict(int))
        for key,val in read_dict_HP1.items():
            hp1_count= val
            hp2_count= read_dict_HP2[key]
            if (hp1_count > hp2_count and 
                hp1_count/(hp1_count+hp2_count) >= hapRatio and 
                hp1_count >= minvariant):
                reads_hap_temp[key]= 1
            elif (hp2_count > hp1_count and 
                  hp2_count/(hp1_count+hp2_count) >= hapRatio and 
                  hp2_count >= minvariant):
                reads_hap_temp[key]= 2
        for key,val in variant_dict.items():
            read_count = hp1_count = hp2_count = hp1_count_read = hp2_count_read = 0
            for read in val:
                read_count += 1
                hp1_count_read += read_dict_HP1[(key[0],*read)]
                hp2_count_read += read_dict_HP2[(key[0],*read)]
                if (key[0],*read) in reads_hap_temp:
                    if reads_hap_temp[(key[0],*read)] == 1:
                        hp1_count += 1
                    elif reads_hap_temp[(key[0],*read)] == 2:
                        hp2_count += 1
            variant_dict_HP1_ave[(key[0],key[1])] = round(hp1_count_read/read_count,3)
            variant_dict_HP2_ave[(key[0],key[1])] = round(hp2_count_read/read_count,3)
            for read in val:
                if (hp1_count > hp2_count and 
                    hp1_count/(hp1_count+hp2_count) >= hapRatio and 
                    hp1_count >= minvariant):
                    read_hap_reassign[(key[0],*read)][1] += 1
                elif (hp2_count > hp1_count and 
                      hp2_count/(hp1_count+hp2_count) >= hapRatio and 
                      hp2_count >= minvariant):
                    read_hap_reassign[(key[0],*read)][2] += 1
        for key in read_hap_reassign.keys():
            hp1_count= read_hap_reassign[key][1]
            hp2_count= read_hap_reassign[key][2]
            if (hp1_count > hp2_count and 
                hp1_count/(hp1_count+hp2_count) >= hapRatio and 
                hp1_count >= minvariant):
                if key[0:2] not in reads_hap:
                    reads_hap[key[0:2]]= 1
                elif key[2] in ["0","16"]:
                    reads_hap[key[0:2]]= 1
            elif (hp2_count > hp1_count and 
                  hp2_count/(hp1_count+hp2_count) >= hapRatio and 
                  hp2_count >= minvariant):
                if key[0:2] not in reads_hap:
                    reads_hap[key[0:2]]= 2
                elif key[2] in ["0","16"]:
                    reads_hap[key[0:2]]= 2
        read_hap_reassign.clear()
        variant_dict_HP1= defaultdict(lambda: defaultdict(int))
        variant_dict_HP2= defaultdict(lambda: defaultdict(int))
        variant_dict_HP1_all= defaultdict(int)
        variant_dict_HP2_all= defaultdict(int)
        for key,val in per_var_cov.items():
            for read in val:
                if read in reads_hap and reads_hap[read]==1:
                    variant_dict_HP1_all[(key[0],key[1])] += 1
                    variant_dict_HP2_all[(key[0],key[1])] += 0
                elif read in reads_hap and reads_hap[read]==2:
                    variant_dict_HP1_all[(key[0],key[1])] += 0
                    variant_dict_HP2_all[(key[0],key[1])] += 1
                
        for key,val in variant_dict.items():
            for read in val:
                if (key[0],read[0]) in reads_hap:
                    if reads_hap[(key[0],read[0])] == 1:
                        variant_dict_HP1_all[(key[0],key[1],"delskip")] += 1
                        variant_dict_HP2_all[(key[0],key[1],"delskip")] += 0
                        variant_dict_HP1[(key[0],key[1])][key[2]] += 1 
                        variant_dict_HP2[(key[0],key[1])][key[2]] += 0
                    elif reads_hap[(key[0],read[0])] == 2:
                        variant_dict_HP1_all[(key[0],key[1],"delskip")] += 0
                        variant_dict_HP2_all[(key[0],key[1],"delskip")] += 1
                        variant_dict_HP1[(key[0],key[1])][key[2]] += 0
                        variant_dict_HP2[(key[0],key[1])][key[2]] += 1 
        variant_dict.clear()
        records_chrom= vcf_tb.query(chrom, 0, chroms[chrom]+1)
        for line in records_chrom:
            block_id= line[9].split(":")[-1]
            if 'PS' in line[8].split(":"):
                ps_index= line[8].split(":").index('PS')
                new_ps= line[8].split(":")
                new_hp= line[9].split(":")
                new_ps.pop(ps_index)
                new_hp.pop(ps_index)
            else:
                new_ps= line[8].split(":")
                new_hp= line[9].split(":")
            
            if not args.include_all_variants and line[6] not in ["PASS","."]:
                continue
            if line[9].startswith(("0/1","1/0","0|1","1|0")):
                hp1_count_alt= variant_dict_HP1[(line[0],int(line[1])-1)][line[4].upper()]
                hp2_count_alt= variant_dict_HP2[(line[0],int(line[1])-1)][line[4].upper()]
                hp1_count_ref= variant_dict_HP1[(line[0],int(line[1])-1)][line[3].upper()]
                hp2_count_ref= variant_dict_HP2[(line[0],int(line[1])-1)][line[3].upper()]
                hp1_cov= str(variant_dict_HP1_all[(line[0],int(line[1])-1)])
                hp2_cov= str(variant_dict_HP2_all[(line[0],int(line[1])-1)])
                hp1_cov_skip= str(variant_dict_HP1_all[(line[0],int(line[1])-1,"delskip")])
                hp2_cov_skip= str(variant_dict_HP2_all[(line[0],int(line[1])-1,"delskip")])
                all_cov= str(len(per_var_cov[(line[0],int(line[1])-1)]))
                if hp1_cov_skip == "0":
                    hp1_count_ave= 0
                else:
                    hp1_count_ave= variant_dict_HP1_ave[(line[0],int(line[1])-1)]
                if hp2_cov_skip == "0":
                    hp2_count_ave= 0
                else:
                    hp2_count_ave= variant_dict_HP2_ave[(line[0],int(line[1])-1)]
                if hp1_count_alt > 0 or hp1_count_ref > 0:
                    hp1_alt_ratio= hp1_count_alt/(hp1_count_alt+hp1_count_ref)
                    hp1_ref_ratio= hp1_count_ref/(hp1_count_alt+hp1_count_ref)
                else:
                    hp1_alt_ratio= 0
                    hp1_ref_ratio= 0
                if hp2_count_alt > 0 or hp2_count_ref > 0:
                    hp2_alt_ratio= hp2_count_alt/(hp2_count_alt+hp2_count_ref)
                    hp2_ref_ratio= hp2_count_ref/(hp2_count_alt+hp2_count_ref)
                else:
                    hp2_alt_ratio= 0
                    hp2_ref_ratio= 0
                if args.phased and (line[0],block_id,"agreement") in phase_block_stat:
                    additional_info= [all_cov,hp1_cov,hp2_cov,hp1_cov_skip,hp2_cov_skip,
                                     str(hp1_count_ref), str(hp2_count_ref),
                                     str(hp1_count_alt),str(hp2_count_alt),
                                     str(hp1_count_ave),
                                     str(hp2_count_ave),
                                     blocks_dict[(line[0],block_id)][0],
                                     blocks_dict[(line[0],block_id)][1],
                                     str(phase_block_stat[(line[0],block_id,"agreement")]),
                                     str(phase_block_stat[(line[0],block_id,"disagreement")])]
                else:
                    additional_info= [all_cov,hp1_cov,hp2_cov,hp1_cov_skip,hp2_cov_skip,
                                     str(hp1_count_ref), str(hp2_count_ref),
                                     str(hp1_count_alt),str(hp2_count_alt),
                                     str(hp1_count_ave),
                                     str(hp2_count_ave),
                                     "NA","NA","NA","NA"]
                if ((hp1_count_alt > hp2_count_alt and
                     hp1_alt_ratio > hp2_alt_ratio and
                     hp1_count_alt >= min_read_reassignment) or
                    (hp2_count_ref > hp1_count_ref and
                     hp2_ref_ratio > hp1_ref_ratio and
                     hp2_count_ref >= min_read_reassignment)):
                    re_assignment_vars[tuple(line[0:2])]= (line[0:8]+
                                                        [':'.join(new_ps)+":PS"]+
                                                        ["1|0:"+':'.join(new_hp[1:])+":HP2|HP1"]+
                                                        additional_info)
                elif ((hp2_count_alt > hp1_count_alt and
                       hp2_alt_ratio > hp1_alt_ratio and
                       hp2_count_alt >= min_read_reassignment) or
                      (hp1_count_ref > hp2_count_ref and
                       hp1_ref_ratio > hp2_ref_ratio and
                       hp1_count_ref >= min_read_reassignment)):
                    re_assignment_vars[tuple(line[0:2])]= (line[0:8]+
                                                          [':'.join(new_ps)+":PS"]+
                                                          ["0|1:"+':'.join(new_hp[1:])+":HP1|HP2"]+
                                                          additional_info)
                else:
                    re_assignment_vars[tuple(line[0:2])]= (line[0:8]+
                                                          [':'.join(new_ps)]+
                                                          [':'.join(new_hp).replace("|", "/").
                                                                            replace("1/0", "0/1")]+
                                                          additional_info)
            elif line[9].startswith(('1/2','1|2','2/1','2|1')):
                hp1_count_alt= variant_dict_HP1[(line[0],int(line[1])-1)][line[4].split(',')[1].upper()]
                hp2_count_alt= variant_dict_HP2[(line[0],int(line[1])-1)][line[4].split(',')[1].upper()]
                hp1_count_ref= variant_dict_HP1[(line[0],int(line[1])-1)][line[4].split(',')[0].upper()]
                hp2_count_ref= variant_dict_HP2[(line[0],int(line[1])-1)][line[4].split(',')[0].upper()]
                hp1_cov= str(variant_dict_HP1_all[(line[0],int(line[1])-1)])
                hp2_cov= str(variant_dict_HP2_all[(line[0],int(line[1])-1)])
                hp1_cov_skip= str(variant_dict_HP1_all[(line[0],int(line[1])-1,"delskip")])
                hp2_cov_skip= str(variant_dict_HP2_all[(line[0],int(line[1])-1,"delskip")])
                all_cov= str(len(per_var_cov[(line[0],int(line[1])-1)]))
                if hp1_cov_skip == "0":
                    hp1_count_ave= 0
                else:
                    hp1_count_ave= variant_dict_HP1_ave[(line[0],int(line[1])-1)]
                if hp2_cov_skip == "0":
                    hp2_count_ave= 0
                else:
                    hp2_count_ave= variant_dict_HP2_ave[(line[0],int(line[1])-1)]
                if hp1_count_alt > 0 or hp1_count_ref > 0:
                    hp1_alt_ratio= hp1_count_alt/(hp1_count_alt+hp1_count_ref)
                    hp1_ref_ratio= hp1_count_ref/(hp1_count_alt+hp1_count_ref)
                else:
                    hp1_alt_ratio= 0
                    hp1_ref_ratio= 0
                if hp2_count_alt > 0 or hp2_count_ref > 0:
                    hp2_alt_ratio= hp2_count_alt/(hp2_count_alt+hp2_count_ref)
                    hp2_ref_ratio= hp2_count_ref/(hp2_count_alt+hp2_count_ref)
                else:
                    hp2_alt_ratio= 0
                    hp2_ref_ratio= 0
                if args.phased and (line[0],block_id,"agreement") in phase_block_stat:
                    additional_info= [all_cov,hp1_cov,hp2_cov,hp1_cov_skip,hp2_cov_skip,
                                     str(hp1_count_ref), str(hp2_count_ref),
                                     str(hp1_count_alt),str(hp2_count_alt),
                                     str(hp1_count_ave),
                                     str(hp2_count_ave),
                                     blocks_dict[(line[0],block_id)][0],
                                     blocks_dict[(line[0],block_id)][1],
                                     str(phase_block_stat[(line[0],block_id,"agreement")]),
                                     str(phase_block_stat[(line[0],block_id,"disagreement")])]
                else:
                    additional_info= [all_cov,hp1_cov,hp2_cov,hp1_cov_skip,hp2_cov_skip,
                                     str(hp1_count_ref), str(hp2_count_ref),
                                     str(hp1_count_alt),str(hp2_count_alt),
                                     str(hp1_count_ave),
                                     str(hp2_count_ave),
                                     "NA","NA","NA","NA"]
                if ((hp1_count_alt > hp2_count_alt and
                     hp1_alt_ratio > hp2_alt_ratio and
                     hp1_count_alt >= min_read_reassignment) or
                    (hp2_count_ref > hp1_count_ref and
                     hp2_ref_ratio > hp1_ref_ratio and
                     hp2_count_ref >= min_read_reassignment)):
                    re_assignment_vars[tuple(line[0:2])]= (line[0:8]+
                                                        [':'.join(new_ps)+":PS"]+
                                                        ["1|2:"+':'.join(new_hp[1:])+":Ref_HP2|HP1"]+
                                                        additional_info)
                        
                elif ((hp2_count_alt > hp1_count_alt and
                       hp2_alt_ratio > hp1_alt_ratio and
                       hp2_count_alt >= min_read_reassignment) or
                      (hp1_count_ref > hp2_count_ref and
                       hp1_ref_ratio > hp2_ref_ratio and
                       hp1_count_ref >= min_read_reassignment)):
                    re_assignment_vars[tuple(line[0:2])]= (line[0:8]+
                                                        [':'.join(new_ps)+":PS"]+
                                                        ["1|2:"+':'.join(new_hp[1:])+":Ref_HP1|HP2"]+
                                                        additional_info)
                else:
                    re_assignment_vars[tuple(line[0:2])]= (line[0:8]+
                                                        [':'.join(new_ps)]+
                                                        [':'.join(new_hp).replace("|", "/").
                                                                          replace("2/1", "1/2")]+
                                                        additional_info)

    subprocess.run("sed '1d' {} | awk -F'\t' '{{print $1,$2-100000,$3+100000}}' OFS='\t'"
                    " | awk '{{if ($2<0) {{$2=0}}; print}}' OFS='\t' "
                    "| sort -k1,1 -k2,2n | bedtools merge > {}"
                    "".format(known_dmr, out+"_temp_knownDMR.tsv"),
                        shell=True,
                        check=True)
    subprocess.run("samtools view -h --remove-tag HP,PS -@ {} "
                   "--regions-file {} {}"
                        " | sed '/^@PG/d' | gzip > {}".format(processes,
                                               out+"_temp_knownDMR.tsv", 
                                               bam_file,
                                                out+"_TempAtDMRs.sam.gz"),
                        shell=True,
                        check=True)
    subprocess.run("rm {}*".format(out+"_temp"),
                    shell=True,
                    check=True)
    
    write_sam_phase(out+"_TempAtDMRs.sam.gz", 
                    out+"_Temp-NonPofO_dmr.sam.gz", 
                    reads_hap)
    subprocess.run("gunzip -c {0} | samtools sort -@ {1} -o {2} && "
                   "samtools index -@ {1} {2}".format(out+"_Temp-NonPofO_dmr.sam.gz",
                                                      processes,
                                                      out+"_Temp-NonPofO_dmr.bam"),
                        shell=True,
                        check=True)
    out_freqhp1= out + '_Temp_NonPofO_HP1-HP2_MethylationHP1.tsv'
    out_freqhp2= out + '_Temp_NonPofO_HP1-HP2_MethylationHP2.tsv'
    out_freq_methbam(out, processes, reference,
                     os.path.abspath(args.pb_cpg_tools_model), 
                     args.pacbio)

    subprocess.run("{} {} {} {} {} {} {} {} {} {} {}".format("Rscript",
                                                os.path.join(os.path.dirname(
                                                            os.path.realpath(__file__)
                                                                ),
                                                          "DMA_UsingDSS.R"),
                                                out_freqhp1,
                                                out_freqhp2,
                                                out+"_Temp",
                                                args.equal_disp,
                                                args.smoothing_flag,
                                                args.smoothing_span,
                                                args.delta_cutoff,
                                                args.pvalue,
                                                dss_processes),
                    shell=True,
                    check=True)
    subprocess.run("bgzip -f {0} && tabix -f -S 1 -p bed {0}.gz"
                    "".format(out+"_Temp_DMLtest.tsv"),
        shell=True,
        check=True)
    subprocess.run("bgzip -f {0} && tabix -f -S 1 -p bed {0}.gz"
                    "".format(out+"_Temp_callDML.tsv"),
        shell=True,
        check=True)
    out_meth= open(out+"_CpG-Methylation-Status-at-DMRs.tsv",'w')
    dmr_file= openfile(known_dmr)
    header= next(dmr_file).rstrip()
    out_meth.write(header+"\tAll_CpGs_At_iDMR_CouldBeExaminedInBothHaplotypes\t"
                   "DifferentiallyMethylatedCpGs_HypermethylatedOnHP1\t"
                   "DifferentiallyMethylatedCpGs_HypermethylatedOnHP2\t"
                   "MethylationFrequency_HP1\tMethylationFrequency_HP2\t"
                   "Included_Or_Ignored_For_PofO_Assignment\n")
    dmr_file.close()
    (chrom_hp_origin_count,
      chrom_hp_origin_count_all_cg,
      chrom_hp_origin_count_diff_cg,
      chrom_hp_origin_count_dmr)= PofO_dmr(known_dmr, 
                                          out,
                                          out_meth,
                                          min_cg,
                                          cpg_difference)
    chrom_hp_origin= pofo_final_dict(chrom_hp_origin_count,
                                     chrom_hp_origin_count_dmr,
                                     chrom_hp_origin_count_all_cg,
                                     chrom_hp_origin_count_diff_cg)


    print("################## Assigning PofO to Variants ##################")
    info_out_dict= defaultdict(lambda: defaultdict(int))
    with openfile(vcf) as vf_file:
        assignment_file= open(out+"_PofO_Assigned.vcf",'w')
        assignment_file_info= open(out+"_PofO_Assigned_info.tsv",'w')
        for line in vf_file:
            if line.startswith("##"):
                assignment_file.write(line)
            elif line.startswith("#"):
                assignment_file.write(line)
                assignment_file_info.write(line.rstrip()+"\tNumAllReads\t"
                                        "NumReadsHP1\tNumReadsHP2\tNumReadsHP1_DelSkipRefSkip\t"
                                        "NumReadsHP2_DelSkipRefSkip\tNumReads_HP1_Ref/Left_Allele"
                                        "\tNumReads_HP2_Ref/Left_Allele\tNumReads_HP1_Alt/Right_Allele"
                                        "\tNumReads_HP2_Alt/Right_Allele\t"
                                        "MeanNumherOfHP1-PhasedVariantsAccrossReads\t"
                                        "MeanNumherOfHP2-PhasedVariantsAccrossReads\t"
                                        "BlockStart\tBlockEnd\t"
                                        "NumberOfSupportiveStrandSeqPhasedVariantsAtThePhasedBlock\t"
                                        "NumberOfConflictingStrandSeqPhasedVariantsAtThePhasedBlock\t"
                                        "NumReads_Maternal_Ref/Left_Allele"
                                        "\tNumReads_Paternal_Ref/Left_Allele\t"
                                        "NumReads_Maternal_Alt/Right_Allele"
                                        "\tNumReads_Paternal_Alt/Right_Allele\t\n")
            else:
                line=line.rstrip().split('\t')
                snv_var= False
                indel_var= False
                if line[9].startswith(("0/1","1/0","0|1","1|0",
                                        "1/2","1|2","2/1","2|1")):
                    if ((len(line[3]) == 1 and len(line[4]) == 1) or 
                        (len(line[3]) == 1 and len(line[4]) == 3 and "," in line[4])):
                        snv_var= True
                        info_out_dict[line[0]]["all_het_snvs"] += 1
                    else:
                        indel_var= True
                        info_out_dict[line[0]]["all_het_indels"] += 1
                if tuple(line[0:2]) in re_assignment_vars:
                    var_info= re_assignment_vars[tuple(line[0:2])]
                    (hp1_count_ref,hp2_count_ref,
                     hp1_count_alt,hp2_count_alt) = var_info[15:19]
                    if line[0] in chrom_hp_origin:
                        if snv_var and var_info[9].startswith(("1|0","0|1","1|2")):
                            info_out_dict[line[0]]["pofo_het_snvs"] += 1
                        elif indel_var and var_info[9].startswith(("1|0","0|1","1|2")):
                            info_out_dict[line[0]]["pofo_het_indels"] += 1
                            
                        if chrom_hp_origin[line[0]]['HP1'][0] == 'maternal':
                            out_line= '\t'.join(var_info[0:10]).replace("HP1", "Mat"
                                                                 ).replace("HP2", "Pat")
                            assignment_file.write(out_line+'\n')
                            assignment_file_info.write(out_line+'\t'+
                                                  '\t'.join(var_info[10:]+
                                                            [hp1_count_ref,hp2_count_ref,
                                                             hp1_count_alt,hp2_count_alt])+'\n')
                        elif chrom_hp_origin[line[0]]['HP1'][0] == 'paternal':
                            if var_info[9].startswith("1|0"):
                                out_line= '\t'.join(var_info[0:10]).replace("1|0","0|1"
                                                                     ).replace("HP1", "Pat"
                                                                     ).replace("HP2", "Mat")
                            elif var_info[9].startswith("0|1"):
                                out_line= '\t'.join(var_info[0:10]).replace("0|1","1|0"
                                                                     ).replace("HP1", "Pat"
                                                                     ).replace("HP2", "Mat")
                            else:
                                out_line= '\t'.join(var_info[0:10]).replace("HP1", "Pat"
                                                                     ).replace("HP2", "Mat")
                            assignment_file.write(out_line+'\n')
                            assignment_file_info.write(out_line+'\t'+
                                                  '\t'.join(var_info[10:]+
                                                            [hp2_count_ref,hp1_count_ref,
                                                             hp2_count_alt,hp1_count_alt])+'\n')
                    else:
                        out_line= '\t'.join(var_info[0:10]).replace(":PS", ""
                                                    ).replace("1|0", "0/1"
                                                    ).replace("2|1", "1/2"
                                                    ).replace("|", "/"
                                                    ).replace(":Ref_HP1/HP2",""
                                                    ).replace(":Ref_HP2/HP1",""
                                                    ).replace(":HP1/HP2",""
                                                    ).replace(":HP2/HP1","")
                        assignment_file.write(out_line+'\n')
                        assignment_file_info.write(out_line+'\t'+
                                              '\t'.join(var_info[10:]+
                                                        ["NA"]*4)+'\n')
                else:
                    if 'PS' in line[8].split(":"):
                        ps_index= line[8].split(":").index('PS')
                        new_ps= line[8].split(":")
                        new_hp= line[9].split(":")
                        new_ps.pop(ps_index)
                        new_hp.pop(ps_index)
                    else:
                        new_ps= line[8].split(":")
                        new_hp= line[9].split(":")
                    out_line= '\t'.join(line[0:8]+[':'.join(new_ps)
                                                   ]+[':'.join(new_hp)]
                                                        ).replace("1|0", "0/1"
                                                        ).replace("2|1", "1/2"
                                                        ).replace("|", "/")
                    assignment_file.write(out_line+'\n')
                    assignment_file_info.write(out_line+'\t'+
                                          '\t'.join(["NA"]*19)+'\n')
    assignment_file.close()
    assignment_file_info.close()

    
    if args.sv_vcf is not None:
        print("################## Assigning PofO to SVs ##################")
        all_sv_files= args.sv_vcf
        for sv_file in all_sv_files:
            pofo_sv_write(sv_file,
                          out,
                          chrom_hp_origin,
                          reads_hap,
                          min_read_reassignment)
    
    print("############ Preparing PofO Tagged Alignment File #############")
    infile_bam = pysam.AlignmentFile(bam_file, "rb")
    cramHeader = infile_bam.header.to_dict()
    if 'PG' in cramHeader:
        cramHeader= cramHeader.pop('PG')
    outfile_bam = pysam.AlignmentFile(out+"_PofO_Tagged.cram", "wc", 
                                      template=infile_bam, header=cramHeader,
                                      reference_filename= reference)
    for chrom in bam_choms:
        alignment_writer(infile_bam,
                         chrom,
                         reads_hap,
                         chrom_hp_origin,
                         outfile_bam)
    outfile_bam.close()
    infile_bam.close()
 
    subprocess.run("rm {}*".format(out+"_Temp"),
                    shell=True,
                    check=True)
                        
    out_scores= open(out + '_PofO_Scores.tsv','w')
    out_scores.write("Chromosome\tOrigin_HP1\tOrigin_HP2\tPofO_Assignment_Score\t"
                     "Num_Differentially_Methylated_CGs_Supported_PofO_Assignment\t"
                     "Num_Differentially_Methylated_CGs_Conflicted_PofO_Assignment\t"
                     "Num_iDMRs_Supported_PofO_Assignment\t"
                     "Num_iDMRs_Conflicted_PofO_Assignment\t"
                     "Num_All_Differentially_Methylated_CGs_At_Supporting_iDMRs\t"
                     "Num_All_Differentially_Methylated_CGs_At_Conflicting_iDMRs\t"
                     "Num_All_CGs_CouldBeExaminedInBothHaplotypes_At_Supporting_iDMRs\t"
                     "Num_All_CGs_CouldBeExaminedInBothHaplotypes_At_Conflicting_iDMRs\n")
    for chrom,val in chrom_hp_origin.items():
        for hp,score in val.items():
            if hp=="HP1":
                origin_hp1= score[0]
                if origin_hp1 == "maternal":
                    origin_hp2= "Paternal"
                elif origin_hp1 == "paternal":
                    origin_hp2= "Maternal"
                out_scores.write('\t'.join([chrom,origin_hp1,origin_hp2,
                                            str(score[1]/(score[1]+score[2]))] +
                                           list(map(str,score[1:])))+'\n')
    out_scores.close()
    if not args.include_all_variants:
        print("Per chromosome info for the variants with \"PASS\""
              " or \".\" in FILTER column in the input vcf file:")
    else:
        print("Per chromosome info for the variants in the input vcf file:")
    print("chrom\tall_het_snvs\tall_het_indels\t"
          "pofo_assigned_het_snvs\tpofo_assigned_het_indels")
    for key,val in info_out_dict.items():
        print('\t'.join(map(str,[key,val["all_het_snvs"],val["all_het_indels"],
                            val["pofo_het_snvs"],val["pofo_het_indels"]])))



"""
Specific argument parser.
"""
parser = argparse.ArgumentParser(prog='patmat.py', add_help=False,
                                 description="Phasing reads and Methylation "
                                             "using strand-seq and nanopore to determine "
                                             "PofO of each homologous chromosome "
                                             "in a single sample.")
required = parser.add_argument_group("Required arguments")
required.add_argument("--bam", "-b",
                      action="store",
                      type=str,
                      required=True,
                      help="The path to the coordinate sorted bam file with methylation"
                      " tag.")
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
                  help="The path to the vcf file. If the input vcf is phased (e.g. using whatshap or longphase) "
                       "you can select the --phased option to also use phase blocks. See --phased for more details.")
required.add_argument("--strand_vcf", "-stv",
                      action="store",
                      type=str,
                      required=True,
                      help="The path to the chromosome-scale phased vcf file."
                           " This is the input vcf file that has been phased "
                           "using strand-seq data.")
required.add_argument("--reference", "-ref",
                      action="store",
                      type=str,
                      required=False,
                      help=("If you have given a bam file with methylation tag"
                            " for the --tool_and_callthresh option, then you "
                            "must also give the path to the reference file. "
                            "File must be indexed using samtools faidx."))
optional = parser.add_argument_group("Optional arguments")
optional.add_argument("--pacbio", "-pb",
                            action="store_true",
                            required=False,
                            help="Select this if the reads are from PacBio HiFi.")

optional.add_argument("--known_dmr", "-kd",
                  action="store",
                  type=str,
                  required= False,
                  default= os.path.join(os.path.dirname(
                                                os.path.realpath(__file__)
                                                    ),
                                             "Imprinted_DMR_List_V1.GRCh38.tsv"),
                  help="The path to the input file for known imprinted DMRs. "
                        "File must have the following information in the "
                        "following column order: "
                        "chromosome\tstart\tend\tMethylatedAlleleOrigin "
                        "where the methylated allele origin must be either "
                        "maternal or paternal (First row must be header). "
                        "By default, we use iDMR list in repo's patmat directory.")
optional.add_argument("--phased", "-ph",
                      action="store_true",
                      required=False,
                      help=("Select this option if your input vcf is a phased vcf file"
                            " using long-reads using whatshap or longphase. "
                            "In this case phased-blocks and phased variants "
                            "inside them will be used and strand-seq phased data"
                            " are used to correct the switches accross phase blocks. "
                            "This can be useful when the chromosome-scale"
                            " phased variants from strand-seq are very sparse. "
                            "If selected, vcf file must be indexed using tabix."
                            " Note that phased blocks will be extract from vcf file and"
                            " assumption is that the number after last \":\" sign in "
                            "column 10 is the block ID, as with vcfs phased by "
                            "whatshap or longphase."))
optional.add_argument("--pb_cpg_tools_model", "-pbcg",
                      action="store",
                            type=str,
                            required=False,
                            default=os.path.join('/'.join(os.path.dirname(
                                                    os.path.realpath(__file__)
                                                        ).split('/')[0:-1]),
                                                 "third_parties",
                                                 "pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu",
                                                 "models/pileup_calling_model.v1.tflite"),
                      help=("If you have selected --pacbio patmat will use "
                            "aligned_bam_to_cpg_scores tool with "
                            "pileup_calling_model.v1.tflite model for methylation"
                            " extraction. You can also use another model using this option."))

optional.add_argument("--sv_vcf", "-sv_vcf",
                      action="store",
                      type=str,
                      nargs='+',
                      required=False,
                      default=None,
                      help=("Path to the Structural variation (SV) call"
                            " file if you wish to also add PofO to SVs."
                            " if multiple files are gived, they must be separated by"
                            " space and give the absolute path to each file."
                            " File(s) must include read names (e.g.RNAMES=read1,"
                            "read2,read3) in the 8th column of vcf ."))
optional.add_argument("--black_list", "-bl",
                      action="store",
                      type=str,
                      required=False,
                      default= None,
                      help="List of regions to ignore phased variants at them."
                      " Three first columns must be chromosome\tstart\tend."
                      " If a black list is given the vcf file must be indexed using tabix.")
optional.add_argument("--hapratio", "-hr",
                      action="store",
                      type=float,
                      required=False,
                      default=0.75,
                      help=("0-1. Minimum ratio of variants a read must have "
                            "from a haplotype to assign it to that haplotype. "
                            "Default is 0.75. Note that if you also provide "
                            "--phased option, this option will be "
                            "also used to correct phased-block switches"
                            " using Strand-seq phased variants. In this case,"
                            " it is the minimum ratio of phased variants "
                            "at a block that supports the decision based "
                            "on strand-seq phased variants."))
optional.add_argument("--mapping_quality", "-mq",
                      action="store",
                      type=int,
                      required=False,
                      default=20,
                      help=("An integer value to specify threshold for "
                            "filtering reads based on mapping quality for "
                            "PofO assignment to variants. "
                            "Default is >=10"))
optional.add_argument("--min_variant", "-mv",
                      action="store",
                      type=int,
                      required=False,
                      default=1,
                      help=("Minimum number of phased variants a read must have "
                            "to be considered during variant rephasing."
                            ". Default= 1."))
optional.add_argument("--min_variant_block", "-mvb",
                      action="store",
                      type=int,
                      required=False,
                      default=2,
                      help=("Minimum number of concordant phased variants must a phased block "
                            "have to be used during phase block switch correction"
                            ". Default= 2."))
optional.add_argument("--min_read_number", "-mr",
                      action="store",
                      type=int,
                      required=False,
                      default=2,
                      help=("Minimum number of reads to support a variant"
                            " to assign to each haplotype. Default= 2"))
optional.add_argument("--min_cg", "-mcg",
                      action="store",
                      type=int,
                      required=False,
                      default= 5,
                      help=("Minimum number of CpGs an iDMR must have to "
                            " consider it for PofO assignment. Default is 5."))
optional.add_argument("--cpg_difference", "-cd",
                      action="store",
                      type=float,
                      required=False,
                      default= 0.1,
                      help=("Minimum cut off for the fraction of CpGs between haplotypes "
                            "must be differentially methylated at an iDMR to "
                            "consider it for PofO assignment. Default is 0.1."))
optional.add_argument("--include_all_variants", "-iav",
                      action="store_true",
                      required=False,
                      help="By default, only variants that have \"PASS\" or \".\" "
                           " in the FILTER column of the input vcf file will be used"
                           " during phasing and PofO assignment. Select this flag "
                           "if you want to use all the variants.")
optional.add_argument("--include_supplementary", "-is",
                      action="store_true",
                      required=False,
                      help="Include supplementary reads.")
optional.add_argument("--processes", "-p",
                      action="store",
                      type=int,
                      required=False,
                      default=4,
                      help="Number of parallel processes. Default is 4.")
optional.add_argument("--chunk_size", "-cs",
                      action="store",
                      type=int,
                      required=False,
                      default=500,
                      help=("Chunk per process. Default is 500"))
optional = parser.add_argument_group("Optional arguments. The following options "
                                     "are DSS options for differential methylation"
                                     " analysis to find differentially methylated "
                                     "CpGs between haplotypes")
optional.add_argument("--delta_cutoff", "-dc",
                            action="store",
                            type=float,
                            default=0.075,
                            required=False,
                            help=("0-1. A threshold for defining differentially "
                                  "methylated loci (DML) or CpGs. In DML testing"
                                  " procedure, hypothesis test that the two groups "
                                  "means are equal is conducted at each CpG site. "
                                  "Here if delta is specified, the function will "
                                  "compute the posterior probability that the "
                                  "difference of the means are greater than delta,"
                                  " and then call DML based on that. Default is 0.075."))
optional.add_argument("--pvalue", "-pv",
                      action="store",
                      type=float,
                      required=False,
                      default= 0.001,
                      help=("0-1. When delta is not specified, this is the "
                            "threshold of p-value for defining DML and "
                            "loci with p-value less than this threshold "
                            "will be deemed DMLs. When delta is specified, "
                            "CpG sites with posterior probability greater than"
                            " 1-pvalue_threshold are deemed DML. Default is 0.001"))
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
                        help=("TRUE/FALSE. A flag to indicate whether to apply "
                              "smoothing in estimating mean methylation levels."
                              " For more instructions see the DSS R package guide. "
                              "Default is TRUE."))
optional.add_argument("--equal_disp", "-ed",
                        action="store",
                        type=str,
                        default="FALSE",
                        required=False,
                        help=("TRUE/FALSE. A flag to indicate whether the "
                              "dispersion in two groups are deemed equal or not. "
                              "For more instructions see the DSS R package guide. "
                              "Default is FALSE Because there is no biological"
                              " replicate here, you should specify either "
                              "equal_disp TRUE or smoothing_flag TRUE. "
                              "Do not specify both as FALSE."))
optional.add_argument("--dss_processes", "-dp",
                      action="store",
                      type=int,
                      required=False,
                      help="Number of parallel processes use for DSS "
                            "differential methylation analysis. If not given, it will "
                            "be the same as --processes option. Differential methylation "
                            " analysis usually requires high memory and if there"
                            " is not enough memory, specify less number processes"
                            " using this flag for DSS to allocate available memory "
                            "for less processes.")
optional = parser.add_argument_group("Help and version options")
optional.add_argument('--version', action='version', 
                      version='%(prog)s 1.4.0_dev',
                      help= "Print program's version and exit")
optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                      help="Print this help and exit.")
args = parser.parse_args()

if __name__ == '__main__':
    main(args)
