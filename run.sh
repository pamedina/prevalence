SRS=${1}
shift 
for SRR in $@; 
do 
#echo "fastq-dump --fasta -I --split-files --stdout ${SRR} >> ${SRS}.fasta 2>> ${SRS}.err"
fastq-dump --fasta -I --split-files --stdout ${SRR} >> ${SRS}.fasta 2>> ${SRS}.err
done 

#echo "./dl_fasta.sh "$*" 2>> ${SRS}.err" 
#./dl_fasta.sh "$*" 2>> ${SRS}.err; 
echo "python3 ./randomly_subsample.new.nopairs.py ${SRS}.fasta > ${SRS}.random.fasta 2>> ${SRS}.err"
python3 ./randomly_subsample.new.nopairs.py ${SRS}.fasta > ${SRS}.random.fasta 2>> ${SRS}.err
echo "rm ${SRS}.fasta"
rm ${SRS}.fasta
echo "blastn -query ${SRS}.random.fasta -db Rickettsiales_4blast.fasta -outfmt 6 > Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast 2>> ${SRS}.err"
blastn -query ${SRS}.random.fasta -db Rickettsiales_4blast.fasta -outfmt 6 > Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast 2>> ${SRS}.err
echo "python3 ./read-rickettsiales-blast.new.py -b Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast -f ${SRS}.random.fasta -e PAIRED > ${SRS}.Rickettsiales_4blast.stats 2>> ${SRS}.err"
python3 ./read-rickettsiales-blast.new.py -b Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast -f ${SRS}.random.fasta -e PAIRED > ${SRS}.Rickettsiales_4blast.stats 2>> ${SRS}.err
#echo "./run_estimate_titer.sh arthropoda_odb9.ancestral.fa-tblastn-${SRS}.blast ${SRS}.random.fasta PAIRED ${SRS} 2>> ${SRS}.err"
#./run_estimate_titer.sh arthropoda_odb9.ancestral.fa-tblastn-${SRS}.random.fasta.blast ${SRS}.random.fasta PAIRED ${SRS} 2>> ${SRS}.err
#echo "./summarize.sh ${SRS} > ${SRS}.arthropoda_odb9.ancestral.summary.stats 2>> ${SRS}.err"
#./summarize.sh ${SRS} > ${SRS}.arthropoda_odb9.ancestral.summary.stats 2>> ${SRS}.err
echo "rm ${SRS}.random.fasta"
rm ${SRS}.random.fasta
