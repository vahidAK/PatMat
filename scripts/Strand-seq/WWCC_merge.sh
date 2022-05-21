for j in "WW" "CC"
do

	for i in $(cut -f1 $1 | sort -u)
	do
		awk -v v1="$i" -v v2="$j" '$1==v1 && $5==v2 {print $2,$3,$4}' OFS="\t" $1 > $i.$j.bed
		samtools view -bh -L $i.$j.bed $i > tmp.$i.$j.bam
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

#could be simplified. Why not take reverse reads and subtract 16, and take forward reads and add 16?
#should also change the mate flags...

samtools view -f64 merged.CC.bam | awk '{$2=$2+64;   print $0}' OFS="\t"  > merged.m1.sam
samtools view -f128 merged.CC.bam | awk '{$2=$2-64;  print $0}' OFS="\t"  > merged.m2.sam

cat <(samtools view -H merged.CC.bam) merged.m1.sam merged.m2.sam |  samtools sort -@$2 -n - | samtools fixmate -@$2 -O BAM - - | samtools sort -@$2 - |
                         samtools view -F4 -F8 -f2 -bh > merged.CC.reversed.bam

else

echo "single end reads"

samtools view -F4 -F8 -F16 merged.CC.bam | awk '{$2=$2+16; $9=-1*$9;  print $0}' OFS="\t"  > merged.f.sam
samtools view -F4 -F8 -f16 merged.CC.bam | awk '{$2=$2-16; $9=-1*$9;  print $0}' OFS="\t"  > merged.r.sam

cat <(samtools view -H merged.CC.bam) merged.f.sam merged.r.sam | samtools sort -@$2 - | samtools view -bh > merged.CC.reversed.bam

fi


rm merged.*.sam

samtools sort -@$2 -n merged.WW.bam | samtools fixmate -O BAM - - | samtools sort -@$2 - |  samtools view -F4 -F8 -f2 -bh| samtools merge -f -@$2 merged.WW.CC.bam merged.CC.reversed.bam -
samtools index -@$2 merged.WW.CC.bam

rm merged.CC.bam* merged.WW.bam* merged.CC.reversed.bam* 
