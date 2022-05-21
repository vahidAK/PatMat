# Edit this header as appropriate, and then run
#	bash /path/to/master.sh	
# From the directory containing your Strand-seq BAM files

#=====================================================================================

# HEADER

# The sample name to use
sample="HG005"

# SNPs called from ONT data
snps="/home/vhanlon/PofO_test/HG005_2-3_Guppy6_Clair3_Passed.vcf"

# Sex of the sample
sex="male"

# The number of threads to use when parallelizing
cpu=24

# Are the reads paired-end? "TRUE" or "FALSE"
paired="TRUE"

# The path to the directory containing these scripts
scripts_dir="/home/vhanlon/PofO_test/scripts/"

# The path to the FASTA-formatted GRCh38 reference genome
reference="/home/vhanlon/sspipe/refseq/human/GRCh38.fasta"

# END OF HEADER

#=====================================================================================

# First we make composite files to call inversions

bash "$scripts_dir"/master_WCCW_composite.sh "$cpu" "$paired" "$scripts_dir" "$snps" "$reference" "$sex" 
bash "$scripts_dir"/master_WWCC_composite.sh "$cpu" "$paired" "$scripts_dir" "$reference" "$sex"

mkdir composite

for j in WW_CC WC_CW
do
	mv "$j"/merged*bam composite/"$sample"."$j".bam
	mv "$j"/merged*bam.bai composite/"$sample"."$j".bam.bai
done

cd composite

Rscript "$scripts_dir"/bpr_simple.R 40 15 "$cpu" "$paired" "$scripts_dir" "BPR_output_small" 
Rscript "$scripts_dir"/bpr_simple.R 120 50 "$cpu" "$paired" "$scripts_dir" "BPR_output_med" 
Rscript "$scripts_dir"/bpr_simple.R 360 50 "$cpu" "$paired" "$scripts_dir" "BPR_output_large"

for i in small med large
do
	cd BPR_output_"$i"/data
	Rscript "$scripts_dir"/extract_bpr.R "$sample" "$i"
	cd ../..
done


cat "$sample"_*bpr.txt | grep -v 'RData' | cut -f1-3 | sort -k1,1 -k2,2n > "$sample".bpr.txt

# Then we actually call the inversions, two different ways, and then combine them

Rscript "$scripts_dir"/runinv.bpr.R "$sample" "$sex" "$paired" "$scripts_dir"
Rscript "$scripts_dir"/runinv.catalogue.R "$sample" "$sex" "$paired" "$scripts_dir"

cat <(bedtools intersect -v -r -f 0.1 -a <(tail -n +2 "$sample".catalogue.genotyped.txt | awk '$9>0.95 && $8!=0 && $8!="0|0"' | sort -k1,1 -k2,2n) -b <(tail -n +2 "$sample".bpr.genotyped.txt | awk '$9>0.95 && $8!=0 && $8!="0|0"' | sort -k1,1 -k2,2n)) <(
		tail -n +2 "$sample".bpr.genotyped.txt | awk '$9>0.95 && $8!=0 && $8!="0|0"' | sort -k1,1 -k2,2n) | sort -k1,1 -k2,2n  > "$sample".inv.txt

awk '$3-$2+1>10000 {print $1,$2,$3}' OFS="\t" "$sample".inv.txt | grep -v '^chrX\|^chrY' > ../"$sample".inv.bed

cd ..

# Next we run BreakpointR and StrandPhaseR to do inversion-aware phasing

Rscript "$scripts_dir"/inversion_phasing.R "$sample" "$snps" "$paired" "$cpu" $(pwd) "$reference" "$scripts_dir"

# Finally we concatenate VCF files from different chromosomes

bcftools concat SPR_output/VCFfiles/chr*INVcorr.vcf > "$sample".phased.inv_aware.vcf
