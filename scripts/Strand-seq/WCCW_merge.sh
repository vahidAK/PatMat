#This script takes the output file from WC_CW_extract.sh as an argument
#It merges the WC regions in all the libraries into a single BAM file; same for the CW regions.
#Afterwards, the CW BAM file has all its reads reversed and is merged with the WC file to create a phase-aware composite file



for j in "WC" "CW"
do

	for i in $(cut -f1 $1 | sort -u)
	do
		awk -v v1="$i" -v v2="$j" '$1==v1 && $5==v2 {print $2,$3,$4}' OFS="\t" $1 > $i.$j.bed
		samtools view -bh -L $i.$j.bed ../$i > tmp.$i.$j.bam
		rm $i.$j.bed
	done

	samtools merge -f -p -c -@$2 merged.$j.bam tmp.*.$j.bam
	samtools index -@$2 merged.$j.bam 
	rm tmp.*.$j.bam

done


#convert R boolean to bash boolean

if [ "$3" = "FALSE" ]
then

	paired=false

elif [ "$3" = "TRUE" ]
then

	paired=true
fi	


#Check whether we're dealing with pe reads here
if $paired
then
	echo "paired end reads"
	#Switching mate one and mate two
	samtools view -f64 merged.CW.bam | awk '{$2=$2+64; print $0}' OFS="\t"  > merged.m1.sam	
	samtools view -f128 merged.CW.bam | awk '{$2=$2-64; print $0}' OFS="\t"  > merged.m2.sam

	cat <(samtools view -H merged.CW.bam) merged.m1.sam merged.m2.sam | samtools sort -@$2 -n - | samtools fixmate -@$2 -O BAM - - | samtools sort -@$2 - |
			samtools view -F4 -F8 -f2 -bh > merged.CW.reversed.bam

else
	echo "single end reads"

	samtools view -F4 -F8 -F16 merged.CW.bam | awk '{$2=$2+16; $9=-1*$9;  print $0}' OFS="\t"  > merged.f.sam
	samtools view -F4 -F8 -f16 merged.CW.bam | awk '{$2=$2-16; $9=-1*$9;  print $0}' OFS="\t"  > merged.r.sam

	cat <(samtools view -H merged.CW.bam) merged.f.sam merged.r.sam | samtools sort -@$2 - | samtools view -bh > merged.CW.reversed.bam


fi


rm merged.*.sam

samtools sort -@$2 -n merged.WC.bam | samtools fixmate -O BAM - - | samtools sort -@$2 - |  samtools view -F4 -F8 -f2 -bh |  samtools merge -f -@$2 merged.WC.CW.bam merged.CW.reversed.bam -
samtools index -@$2 merged.WC.CW.bam

mv merged.WC.CW.bam* ..

rm merged.WC.bam* merged.CW.bam* merged.CW.reversed.bam* 

