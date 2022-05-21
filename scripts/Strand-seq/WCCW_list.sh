#Run this in a StrandPhaseR output directory and it will distinguish WC from CW segments by library.
#The Phased/ directory contains tables exactly this information, so this script is just processing text files.



>wc_cw_strand_state.txt


cd Phased

for i in chr*_phasedFiles_hap1.txt
do
	chr=$(echo $i | cut -f1 -d_)
	tail -n +2 $i > $i.tmp 

	paste <(cut -f1 -d: $i.tmp | rev | cut -f3- -d_ | rev | sed 's/\"//g')  <(cut -f2 -d: $i.tmp | cut -f1 -d-) <(cut -f2 -d: $i.tmp | cut -f2 -d- |  cut -f1 -d_) \
	<(sed 's/_C\"/_CW/g' $i.tmp | sed 's/_W\"/_WC/g' | cut -f1 -d" " | rev | cut -f1 -d_  | rev) \
	| awk -v var="$chr" '{print $1, var, $2, $3, $4}' OFS="\t"  >> ../wc_cw_strand_state.txt
	
	rm $i.tmp

done

cd ..
