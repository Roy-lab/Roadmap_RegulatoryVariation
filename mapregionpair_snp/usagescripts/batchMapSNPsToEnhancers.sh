INDIR=../../featurecombinations/10times/top20/outputs/genomewidepred

for CELL in  H1hesc Helas K562 Gm12878 Hepg2
#for CELL in Gm12878 
do
	#for FILE in Maurano GWAS
	for FILE in GWAS #eqtl
	do
	#./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/../../../../genomewide/${CELL}_pair_nw.txt $INDIR/$CELL/all_enhancer_genematch_${FILE} E
	#./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/../../../../genomewide/${CELL}_pair_promoter_nw.txt $INDIR/$CELL/all_promoter_genematch_${FILE} E
	#./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/$CELL/merged_percentile_avg_top1_nw.txt $INDIR/$CELL/top1_enhancer_genematch_${FILE} E
	./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/$CELL/merged_percentile_avg_top1_nw.txt $INDIR/$CELL/top1_enhancer_genematch_${FILE} E
	HIT=`cat $INDIR/${CELL}/top1_enhancer_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}'  | egrep -v REGULATORY| awk '{printf("%s\t%s\n",$3,$7)}' | sort -u |awk '$2>0'  |wc -l`
	TOTAL=`cat $INDIR/${CELL}/top1_enhancer_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}' | egrep -v REGULATORY |awk '{printf("%s\t%s\n",$3,$7)}' | sort -u |awk '$2>=0'  |wc -l`
	REGULATORY=`cat $INDIR/${CELL}/top1_enhancer_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}'  | egrep REGULATORY | awk '{printf("%s\t%s\n",$3,$7)}' | sort -u   |wc -l`
	#echo "$INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt"
	CODING=`cat $INDIR/${CELL}/top1_enhancer_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}' | egrep  CODING |awk '{printf("%s\t%s\n",$3,$7)}' | sort -u  |wc -l`
	echo  "$CELL $REGULATORY $CODING $HIT $TOTAL" | awk '{printf("%s\t%d\t%d\t%d\t%d\t%f\n",$1,($2+$3),$2,$3,$4,$4/$5)}'
	#echo  "$CELL $HIT $TOTAL" | awk '{printf("%s\t%d\t%d\t%f\n",$1,$2,$3,($2/$3))}'
	#./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/$CELL/merged_percentile_avg_top1_promoter_nw.txt $INDIR/$CELL/top1_promoter_genematch_${FILE} E
	#./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/$CELL/merged_percentile_avg_top1_nw.txt tmp
	done
done
