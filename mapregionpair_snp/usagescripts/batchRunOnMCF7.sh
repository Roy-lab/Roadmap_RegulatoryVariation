DATADIR=/mnt/ws/sysbio/roygroup/shared/projects/e_p_project/mcf7_highthroughput_hic/formatted
INPUT=../../../breastcancer/data/KnowEng/ERpos_LDSnps/
OUTPUT=../../../breastcancer/results/
for D in 5C_MCF7 gisChiaPet_MCF7_Ctcf gisChiaPet_MCF7_Eraa gisChiaPet_MCF7_Pol2 gisChiaPet_MCF7_RNAPII_pilot gisChiaPet_MCF7_RNAPII_saturated 
do
	for LD in 0.1 0.8
	do
	for Rtype in E P
	do
		INFNAME=$INPUT/ERpos_LD_SNPs_r2_${LD}_no_indels_reformat.txt
		PAIRFNAME=${DATADIR}/${D}.txt
		./mapEnhancerToGwas ${INFNAME} null $PAIRFNAME  $OUTPUT/ERpos_SNPs_${LD}_${D}_${Rtype} 	${Rtype}
		echo ">>>>Result of Snp intersect with $D $Rtype"
	done
	done
done
