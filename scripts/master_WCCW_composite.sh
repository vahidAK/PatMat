#Master script that creates a phased WCCW BAM file from a set of Strand-seq libraries for a sample


############################## User sets variables here ###################################################


#parallelize the composite file creation
threads="$1"

#are the reads in the BAM files paired-end?
paired="$2" 

#path to the directory containing these scripts
scripts="$3"

#The path to a VCF file containing good-quality heterozygous SNPs for the sample. If empty, this script will try to call SNPs with bbtools.
snps="$4"

#reference genome
ref="$5"

#sex of the sample
sex="$6"


##########################################################################################################
##########################################################################################################
################################# StrandPhaseR ###########################################################

#WC/CW regions need to be treated separately, which is a multi-step process
#StrandPhaseR will be used, and that requires BreakpointR output (run here with windowsize=6000, minRead=50)
Rscript $scripts/bpr.R 20000000 50 $threads $paired $sex "WC" $scripts TRUE

mkdir ./WC_CW
Rscript $scripts/printWCCWregions.R

cd ./WC_CW
sed 's/:/ /g' wc_regions.txt | awk '{print $1,$2,$3,$4}' OFS="\t" > tmp
mv tmp wc_regions.txt


if [ "$snps" = "" ]; then

	bash $scripts/call_SNPs.sh $ref	
	snps=$(pwd)"/called_snps.vcf"
fi

bcftools view -H -v snps -m2 -M2 -g het $snps > snps.vcf
snps=$(pwd)"/snps.vcf"

#Calling StrandPhaseR
Rscript $scripts/strandPhase.R $threads $(pwd) $snps $paired

################################ Custom WC CW composite file ###############################################

cd SPR_output

#StrandPhaseR lists WC/CW regions by BAM file in the Phased/ output directory, and these can be extracted
bash $scripts/WCCW_list.sh

cd ..
mv SPR_output/wc_cw_strand_state.txt .

#Merging WC and CW regions of BAM files, then reversing all reads in the CW BAM file to merge it with the WC reads
bash $scripts/WCCW_merge.sh wc_cw_strand_state.txt $threads $paired


#An optional step to manually check the file ---- BPR results are useful for that. 
mv ../merged.WC.CW.bam* .
Rscript $scripts/bpr.R 100000 50 $threads $paired $sex "WC" $scripts FALSE
