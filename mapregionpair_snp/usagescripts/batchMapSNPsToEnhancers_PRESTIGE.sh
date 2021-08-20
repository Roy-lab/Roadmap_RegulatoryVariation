INDIR=../../data/prestige

for CELL in  H1hesc Helas K562 Gm12878 Hepg2
do
	#for FILE in eqtl #GWAS
	for FILE in GWAS
	do
	./mapEnhancerToGwas ~/data/human/gwas/${FILE}_formated.txt ~/data/human/genenames/hgnc_complete_set.txt $INDIR/${CELL}_High_enhancer_gene.txt ${INDIR}/${CELL}_genematch_${FILE} E
	HIT=`cat $INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}'  | egrep -v REGULATORY| awk '{printf("%s\t%s\n",$3,$7)}' | sort -u |awk '$2>0'  |wc -l`
	#echo "$INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt"
	TOTAL=`cat $INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}' | egrep -v REGULATORY |awk '{printf("%s\t%s\n",$3,$7)}' | sort -u |awk '$2>=0'  |wc -l`
	#echo  "$CELL $HIT $TOTAL" | awk '{printf("%s\t%d\t%d\t%f\n",$1,$2,$3,$2/$3)}'
	REGULATORY=`cat $INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}'  | egrep REGULATORY | awk '{printf("%s\t%s\n",$3,$7)}' | sort -u   |wc -l`
	#echo "$INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt"
	CODING=`cat $INDIR/${CELL}_genematch_${FILE}_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}' | egrep  CODING |awk '{printf("%s\t%s\n",$3,$7)}' | sort -u  |wc -l`
	echo  "$CELL $REGULATORY $CODING $HIT $TOTAL" | awk '{printf("%s\t%d\t%d\t%d\t%d\t%f\n",$1,($2+$3),$2,$3,$4,$4/$5)}'
	done
done
