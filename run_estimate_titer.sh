BLAST=$1
FASTA=$2
END=$3
SRR=$4
for escore in 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1; 
do 
python3 ./estimate_titer.final.py -b $BLAST -f $FASTA -e $END -m $escore > arthropoda_odb9.ancestral.fa-tblastn-"$SRR".escore"$escore".out 2> arthropoda_odb9.ancestral.fa-tblastn-"$SRR".escore"$escore".out.err
done
