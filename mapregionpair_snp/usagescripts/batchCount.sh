for CELL in H1hesc Helas K562 Gm12878 Hepg2
#for CELL in H1hesc K562 Gm12878 Hepg2
do
echo "$CELL"
HIT=`cat ../../data/prestige/${CELL}_genematch_GWAS_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}'  | egrep -v REGULATORY| awk '{printf("%s\t%s\n",$3,$7)}' | sort -u |awk '$2>0'  |wc -l`
TOTAL=`cat ../../data/prestige/${CELL}_genematch_GWAS_enhancer_snp.txt  | awk '{gsub(/ /,"_"); print $0}' | egrep -v REGULATORY |awk '{printf("%s\t%s\n",$3,$7)}' | sort -u |awk '$2>=0'  |wc -l`
echo  "$CELL $HIT $TOTAL"
done
