OUTDIR=/work/sroy/enhancer_promoter/featurecombinations/distcontrolled_v2/rf_v2_noK27_p300/genomewidepred/universal_enhancer_promoter
OUTDIR=/work/sroy/enhancer_promoter/featurecombinations/distcontrolled_v2/rf_v3_glasso_intersect_nok4me3_k79me2_Smc3_yesk20me1_Tbp/
for CELL in K562 Gm12878 Hepg2 H1hesc Helas
do
	if [ $CELL == Helas ]
	then
		PCELL=Helas3
	else
		PCELL=$CELL
	fi
	#SUFF=percentile_sorted_0.9
	SUFF=percentile_random_0.9
	echo "./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_sorted_top1_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt proximal $OUTDIR/$CELL/rfNames_${SUFF}_sorted_top1_proximal_nw.txt"
	./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt proximal $OUTDIR/$CELL/rfNames_${SUFF}_proximal_nw_deg1.txt
	echo "./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_sorted_top1_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt distal $OUTDIR/$CELL/rfNames_${SUFF}_sorted_top1_distal_nw.txt"
	./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt distal $OUTDIR/$CELL/rfNames_${SUFF}_distal_nw_deg1.txt
	#echo "./genDistalProxNet $OUTDIR/$CELL/merged_percentile_avg_bottom1_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt distal $OUTDIR/$CELL/bottom1_distal_nw.txt"
	#./genDistalProxNet $OUTDIR/$CELL/merged_percentile_avg_bottom1_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt distal $OUTDIR/$CELL/bottom1_distal_nw.txt
done
