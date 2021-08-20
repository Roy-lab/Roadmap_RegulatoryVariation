OUTDIR=/work/sroy/enhancer_promoter/featurecombinations/distcontrolled_v2/rf_v2_noK27_p300/genomewidepred/universal_enhancer_promoter
OUTDIR=../..//featurecombinations/distcontrolled_v2/1times/set1/multicellline_sig_rf_v2_noK27_p300_features_bin_cc_proper/train/outputs/
#for CELL in K562 Gm12878 Hepg2 H1hesc Helas
for CELL in K562 Gm12878 H1hesc Helas
do
	if [ $CELL == Helas ]
	then
		PCELL=Helas3
	else
		PCELL=$CELL
	fi
	for SUFF in enhanceronly enhanceronly_pred0.5
	do
		echo "./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt proximal $OUTDIR/$CELL/rfNames_${SUFF}_proximal_nw.txt"
		./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt proximal $OUTDIR/$CELL/rfNames_${SUFF}_proximal_nw.txt
		echo "./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt distal $OUTDIR/$CELL/rfNames_${SUFF}_distal_nw.txt"
		./genDistalProxNet $OUTDIR/$CELL/rfNames_${SUFF}_genenames.txt /mnt/ws/sysbio/roygroup/shared/data_new/human/encode_narrowpeak/dnase_narrowpeak/idrPool.wgEncodeOpenChromDnase${PCELL}Rep0.narrowPeak ../../data/regnet/dnase1_motifs/wgEncodeOpenChromDnase${PCELL}Rep0_FIMO_all.txt distal $OUTDIR/$CELL/rfNames_${SUFF}_distal_nw.txt
	done
done
