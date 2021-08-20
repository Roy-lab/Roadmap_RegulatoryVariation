DATADIR=../../data/chia-pet/snyder/
for D in GM12878_RAD21.txt K562_H4K27ac.txt K562_H4K4Me1.txt K562_H4K4Me2.txt K562_H4K4Me3.txt K562_PolII.txt K562_RAD21.txt
do
	for Rtype in E P
	do
	./mapEnhancerToGwas ../../data/KnowEng/SNP_Region0_reformat.txt ~/data/human/genenames/hgnc_complete_set.txt ../../data/chia-pet/snyder/K562_H4K27ac.txt Knoweng_temp ${Rtype}	> /dev/null
	echo ">>>>Result of Snp intersect with $D $Rtype"
	cat Knoweng_temp_enhancer_snp.txt
	done
done
