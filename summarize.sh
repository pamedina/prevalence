SRR=$1
for ESCORE in 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1;
do
#echo "$ESCORE" 
TOTAL_BLAST_HITS=$(grep -v nan arthropoda_odb9.ancestral.fa-tblastn-"$SRR".escore"$ESCORE".out | awk '{SUM += $5} END {print SUM}')
GENE_LENGTH=$(grep -v nan arthropoda_odb9.ancestral.fa-tblastn-"$SRR".escore"$ESCORE".out | awk '{SUM += $6} END {print SUM}')
TOTAL_NUMBER_OF_READS_HIT=$(grep -v nan arthropoda_odb9.ancestral.fa-tblastn-"$SRR".escore"$ESCORE".out | awk '{SUM += $13} END {print SUM}')
#COVERAGE=$(($TOTAL_BLAST_HITS/$GENE_LENGTH))
#echo TOTAL_BLAST_HITS $TOTAL_BLAST_HITS
#echo GENE_LENGTH $GENE_LENGTH

PROP=$(bc -l <<< "$TOTAL_BLAST_HITS/$GENE_LENGTH")
#echo $COVERAGE
AVERAGE_OF_MAX_DEPTH_ACROSS_GENES=$(grep -v nan arthropoda_odb9.ancestral.fa-tblastn-"$SRR".escore"$ESCORE".out | awk '{SUM += $12} END {print SUM/NR}')
echo tblastn "$SRR" "$ESCORE" "$TOTAL_BLAST_HITS" "$GENE_LENGTH" "$PROP" "$AVERAGE_OF_MAX_DEPTH_ACROSS_GENES" "$TOTAL_NUMBER_OF_READS_HIT" 
done
