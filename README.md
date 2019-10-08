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

    $ mkdir -p ~/.ncbi
    $ echo '/repository/user/main/public/root = "/scratch/standage/sra-cache"' > ~/.ncbi/user-settings.mkfg
    $ echo '/repository/user/cache-disabled = "true"' > ~/.ncbi/user-settings.mkfg
    
#### Note: 

Information on our method's accuracy to identify divergent strains can be found in Medina et al. (2019) on the bioRxiv. 

Ther version used to produce the results in Medina et al. (2019) is maintained as version_0.1/

#### Basic Usage: 
    $ ./run.sh {SRS} {SRR1} {SRR2} {SRR3} ... 

#### Basic Workflow: 
    $ fastq-dump --fasta -I --split-files --stdout ${SRR} >> ${SRS}.fasta 2>> ${SRS}.err  
    $ python3 randomly_subsample.new.nopairs.py ${SRS}.fasta > ${SRS}.random.fasta 2>> ${SRS}.err  
    $ blastn -query ${SRS}.random.fasta -db Rickettsiales_4blast.fasta -outfmt 6 > Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast 2>> ${SRS}.err  
    $ python3 ./read-rickettsiales-blast.new.py -b Rickettsiales_4blast.fasta-blastn-${SRS}.random.fasta.blast -f ${SRS}.random.fasta -e PAIRED > ${SRS}.Rickettsiales_4blast.stats 2>> ${SRS}.err  

#### Output file format 
The output file from `read-rickettsiales-blast.new.py` is tab delimited. Each sample will have a separate out file, taking the form of ${SRS}.Rickettsiales_4blast.stats, as outlined above. The columns of the file are described below:  
	1. Sample name.  
	2. Human readable reference name.  
	3. Reference name.  
	4. End spots.  
	5. Download number.  
	6. Number of blast hits (bp).  
	7. Length of reference (bp).  
	8. Estimated coverage of sample.  
	9. Variance coefficient of depth.  
	10. Proportion of bins filled.  
	11. Raw bin counts.  
	12. Number of reads mapping to reference.  


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

