DATADIR=../../data/chia-pet/
for D in  K562_RNAPII.txt K562_CTCF.txt
do
	for Rtype in E P
	do
	./mapEnhancerToGwas ../../data/KnowEng/SNP_Region0_reformat.txt ~/data/human/genenames/hgnc_complete_set.txt ../../data/chia-pet/snyder/K562_H4K27ac.txt Knoweng_temp ${Rtype}	> /dev/null
	echo ">>>>Result of Snp intersect with $D $Rtype"
	cat Knoweng_temp_enhancer_snp.txt
	done
done
