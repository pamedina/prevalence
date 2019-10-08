# MICROBE SCAN 

### Paloma  Medina 

#### Download and Compile: 

    $ git clone XXX
    $ cd microbe_scan/src/
    $ make 
  
#### Dependencies 

The software requires blastn and sra-toolkit.  
More information and detailed instructions for blast can be found [here.](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
And information to download the SRA toolkit can be found [here.](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)

    $ sudo apt-get install ncbi-blast+  
    $ sudo apt-get install sra-toolkit  
  
It is also recommended that users disable the SRA cache system as downloading sequence data from the SRA.   
Disabling the SRA cache system will save disk space.   

#### Note: 

Information on our method's accuracy to identify divergent strains can be found in Medina et al. (2019) on the bioRxiv. 

Ther version used to produce the results in Medina et al. (2019) is maintained as version_0.1/

#### Basic Usage:
  $ fastq-dump --fasta -I --split-files --stdout ${SRR} >> ${SRS}.fasta 2>> ${SRS}.err  
  $ python3 randomly_subsample.new.nopairs.py ${SRS}.fasta > ${SRS}.random.fasta 2>> ${SRS}.err  
  $ blastn -query ${SRS}.random.fasta -db Rickettsiales_4blast.fasta -outfmt 6 > Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast 2>> ${SRS}.err  
  $ python3 ./read-rickettsiales-blast.new.py -b Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast -f ${SRS}.random.fasta -e PAIRED > ${SRS}.Rickettsiales_4blast.stats 2>> ${SRS}.err  
  
#### Input file format 


#### Output file format 


#### Original Article 
Medina, P., Russell, S.L. and Corbett-Detig, R., 2019. Deep data mining reveals variable abundance and distribution of microbial reproductive manipulators within and among diverse host species. bioRxiv, https://doi.org/10.1101/679837.

@article {Medina679837,
	author = {Medina, Paloma and Russell, Shelbi L and Corbett-Detig, Russell},
	title = {Deep data mining reveals variable abundance and distribution of microbial reproductive manipulators within and among diverse host species},
	elocation-id = {679837},
	year = {2019},
	doi = {10.1101/679837},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2019/06/23/679837},
	eprint = {https://www.biorxiv.org/content/early/2019/06/23/679837.full.pdf},
	journal = {bioRxiv}
}


# prevalence
scripts associated with "Deep data mining reveals variable abundance and distribution of microbial reproductive manipulators within and among diverse host species"

