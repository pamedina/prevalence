#!/usr/bin/env/ python3 
#Paloma MEdina (pamedina)
#2018-08-29
#group: none 

"""
This script summarizes a blast file. 
The blastdatabase file is /public/home/pamedina/corbettlab/pamedina/WOLBACHIA/AXEL/REFERENCE/Wol_genomes_4blast.fasta.cmb

command line usage : 
python3 read-rickettsiales-blast.py -b SRR035504.fasta.blast -f SRR035504.fasta -s 2437477 -e pair
"""

import sys
import os
import datetime 

WOL_GENOME_SIZES={"NZ_CP029619.1":1193042, "GCF_003788695.1_ASM378869v1_genomic.merged.fna":1358212, "NC_018605.1":944930, "NZ_CBQZ01000000X":1012588, "NZ_CP022339.1":1103593,"NC_002978.6": 1267782,"NC_006833.1":1080084,"NZ_DS996944.1":1543661,"NC_018267.1":957990,"NZ_AP013028.1":1250060,"NZ_MJMG01000001.1":975127,"GCF_000007025.1_ASM702v1_genomic.fna.cmb":1268755,"GCF_000008045.1_ASM804v1_genomic.fna.cmb":1111496,"GCF_000012145.1_ASM1214v1_genomic.fna.cmb":1587240,"GCF_000012385.1_ASM1238v1_genomic.fna.cmb":1522076,"GCF_000014345.1_ASM1434v1_genomic.fna.cmb":1159772,"GCF_000016625.1_ASM1662v1_genomic.fna.cmb":1376184,"GCF_000017445.4_ASM1744v3_genomic.fna.cmb":1268201,"GCF_000018205.1_ASM1820v1_genomic.fna.cmb":1231060,"GCF_000018225.1_ASM1822v1_genomic.fna.cmb":1257710,"GCF_000018245.1_ASM1824v1_genomic.fna.cmb":1528980,"GCF_000021525.1_ASM2152v1_genomic.fna.cmb":1314898,"GCF_000022785.1_ASM2278v1_genomic.fna.cmb":1111612,"GCF_000023005.1_ASM2300v1_genomic.fna.cmb":1290917,"GCF_000195735.1_ASM19573v1_genomic.fna.cmb":1111523,"GCF_000221205.1_ASM22120v1_genomic.fna.cmb":1278471,"GCF_000237845.1_ASM23784v1_genomic.fna.cmb":1275089,"GCF_000252365.1_ASM25236v1_genomic.fna.cmb":1275720,"GCF_000277165.1_ASM27716v1_genomic.fna.cmb":1109804,"GCF_000277185.1_ASM27718v1_genomic.fna.cmb":1111454,"GCF_000277205.1_ASM27720v1_genomic.fna.cmb":1111445,"GCF_000277245.1_ASM27724v1_genomic.fna.cmb":1111969,"GCF_000277265.1_ASM27726v1_genomic.fna.cmb":1112101,"GCF_000277285.1_ASM27728v1_genomic.fna.cmb":1112372,"GCF_000277305.1_ASM27730v1_genomic.fna.cmb":1112957,"GCF_000283595.1_ASM28359v1_genomic.fna.cmb":1283087,"GCF_000283775.1_ASM28377v1_genomic.fna.cmb":1270083,"GCF_000283795.1_ASM28379v1_genomic.fna.cmb":1267197,"GCF_000283815.1_ASM28381v1_genomic.fna.cmb":1269837,"GCF_000283915.1_ASM28391v1_genomic.fna.cmb":1150228,"GCF_000283935.1_ASM28393v1_genomic.fna.cmb":1270751,"GCF_000283955.1_ASM28395v1_genomic.fna.cmb":1255681,"GCF_000283995.1_ASM28399v1_genomic.fna.cmb":1287740,"GCF_000284055.1_ASM28405v1_genomic.fna.cmb":1480884,"GCF_000284075.1_ASM28407v1_genomic.fna.cmb":1305467,"GCF_000284155.1_ASM28415v1_genomic.fna.cmb":1323280,"GCF_000284175.1_ASM28417v1_genomic.fna.cmb":1279798,"GCF_000284195.1_ASM28419v1_genomic.fna.cmb":1300386,"GCF_000363905.1_ASM36390v1_genomic.fna.cmb":1111520,"GCF_000367405.1_ASM36740v1_genomic.fna.cmb":1109301,"GCF_000499665.2_RMONA_1_genomic.fna.cmb":1353450,"GCF_000831525.1_ASM83152v1_genomic.fna.cmb":1257005,"GCF_000831545.1_ASM83154v1_genomic.fna.cmb":1269809,"GCF_001273795.1_ASM127379v1_genomic.fna.cmb":1456213,"GCF_001442475.1_ASM144247v1_genomic.fna.cmb":1448632,"GCF_001602215.1_ASM160221v1_genomic.fna.cmb":1111769,"GCF_001950995.1_ASM195099v1_genomic.fna.cmb":1268220,"GCF_001951015.1_ASM195101v1_genomic.fna.cmb":1268242,"GCF_002078315.1_ASM207831v1_genomic.fna.cmb":1579798,"GCF_002078335.1_ASM207833v1_genomic.fna.cmb":1480876,"GCF_002285905.1_ASM228590v1_genomic.fna.cmb":1378618,"GCF_002356695.1_ASM235669v1_genomic.fna.cmb":1284030,"GCF_002357115.1_ASM235711v1_genomic.fna.cmb":1284030,"GCF_002357135.1_ASM235713v1_genomic.fna.cmb":1283942,"GCF_002357155.1_ASM235715v1_genomic.fna.cmb":1283227,"GCF_002357175.1_ASM235717v1_genomic.fna.cmb":1284030,"GCF_002357195.1_ASM235719v1_genomic.fna.cmb":1284037,"GCF_002357215.1_ASM235721v1_genomic.fna.cmb":1284030,"GCF_002357235.1_ASM235723v1_genomic.fna.cmb":1284009,"GCF_002357255.1_ASM235725v1_genomic.fna.cmb":1284028,"GCF_002357275.1_ASM235727v1_genomic.fna.cmb":1283942,"GCF_002357295.1_ASM235729v1_genomic.fna.cmb":1284030,"GCF_000429565.1_ASM42956v1_genomic.fna.cmb":3670548,"GCF_000757905.1_ASM75790v1_genomic.fna.cmb":2953863,"GCF_001534665.1_ASM153466v1_genomic.fna.cmb":836724,"GCF_000400935.1_ASM40093v1_genomic.fna.cmb":1123322,"GCF_000400955.1_ASM40095v1_genomic.fna.cmb":1107344,"GCF_000439435.1_ASM43943v1_genomic.fna.cmb":1086278,"GCF_000439455.1_ASM43945v1_genomic.fna.cmb":945296,"GCF_000500935.1_ASM50093v1_genomic.fna.cmb":1160554,"GCF_000517365.1_ASM51736v1_genomic.fna.cmb":1132591,"GCF_000565175.1_ASM56517v1_genomic.fna.cmb":1175131,"GCF_000565195.1_ASM56519v1_genomic.fna.cmb":1132608,"GCF_000565215.1_ASM56521v1_genomic.fna.cmb":1075953,"GCF_001029245.1_ASM102924v1_genomic.fna.cmb":1160484,"GCF_001029265.1_ASM102926v1_genomic.fna.cmb":1365714,"GCF_001262715.1_ASM126271v1_genomic.fna.cmb":1261374,"GCF_001267155.1_ASM126715v1_genomic.fna.cmb":1225519,"GCF_001274875.1_ASM127487v1_genomic.fna.cmb":1514561,"GCF_001281045.1_ASM128104v1_genomic.fna.cmb":1179577,"GCF_001513715.1_ASM151371v1_genomic.fna.cmb":1261375,"GCF_001715535.1_ASM171553v1_genomic.fna.cmb":1326546,"GCF_001792795.1_ASM179279v1_genomic.fna.cmb":1199640,"GCF_001886495.1_ASM188649v1_genomic.fna.cmb":1199621,"GCF_001886855.1_ASM188685v1_genomic.fna.cmb":1640878,"GCF_002028345.1_ASM202834v1_genomic.fna.cmb":1364757,"GCF_002237575.1_ASM223757v1_genomic.fna.cmb":1204639,"GCF_002795265.1_ASM279526v1_genomic.fna.cmb":1560007,"GCF_002813555.1_ASM281355v1_genomic.fna.cmb":1284130,"GCF_002865545.1_ASM286554v1_genomic.fna.cmb":891575,"GCF_003339775.1_ASM333977v1_genomic.fna.cmb":1908276,"GCF_003363775.1_ASM336377v1_genomic.fna.cmb":1336077,"GCF_000953015.1_Candidatus_Methylopumilus_turicensis_MMS-10A-171_genomic.fna.cmb":1754988,"GCF_000981505.1_Candidatus_Methylopumilus_planktonicus_MMS-2-53_genomic.fna.cmb":1356428}

COMMON_NAME = {"NZ_CP029619.1":"cardinium1", "GCF_003788695.1_ASM378869v1_genomic.merged.fna":"cardinium2", "NC_018605.1":"cardinium3", "NZ_CBQZ01000000X":"cardinium4", "NZ_CP022339.1":"cardinium5","GCF_000429565.1_ASM42956v1_genomic.fna.cmb":"arsenophonus1","GCF_001534665.1_ASM153466v1_genomic.fna.cmb":"arsenophonus2","GCF_000757905.1_ASM75790v1_genomic.fna.cmb":"arsenophonus3","GCF_000400935.1_ASM40093v1_genomic.fna.cmb":"spriplasma1","GCF_001281045.1_ASM128104v1_genomic.fna.cmb":"spriplasma2","GCF_000400955.1_ASM40095v1_genomic.fna.cmb":"spriplasma3","GCF_001513715.1_ASM151371v1_genomic.fna.cmb":"spriplasma4","GCF_000439435.1_ASM43943v1_genomic.fna.cmb":"spriplasma5","GCF_001715535.1_ASM171553v1_genomic.fna.cmb":"spriplasma6","GCF_000439455.1_ASM43945v1_genomic.fna.cmb":"spriplasma7","GCF_001792795.1_ASM179279v1_genomic.fna.cmb":"spriplasma8","GCF_000500935.1_ASM50093v1_genomic.fna.cmb":"spriplasma9","GCF_001886495.1_ASM188649v1_genomic.fna.cmb":"spriplasma10","GCF_000517365.1_ASM51736v1_genomic.fna.cmb":"spriplasma11","GCF_001886855.1_ASM188685v1_genomic.fna.cmb":"spriplasma12","GCF_000565175.1_ASM56517v1_genomic.fna.cmb":"spriplasma13","GCF_002028345.1_ASM202834v1_genomic.fna.cmb":"spriplasma14","GCF_000565195.1_ASM56519v1_genomic.fna.cmb":"spriplasma15","GCF_002237575.1_ASM223757v1_genomic.fna.cmb":"spriplasma16","GCF_000565215.1_ASM56521v1_genomic.fna.cmb":"spriplasma17","GCF_002795265.1_ASM279526v1_genomic.fna.cmb":"spriplasma18","GCF_001029245.1_ASM102924v1_genomic.fna.cmb":"spriplasma19","GCF_002813555.1_ASM281355v1_genomic.fna.cmb":"spriplasma20","GCF_001029265.1_ASM102926v1_genomic.fna.cmb":"spriplasma21","GCF_002865545.1_ASM286554v1_genomic.fna.cmb":"spriplasma22","GCF_001262715.1_ASM126271v1_genomic.fna.cmb":"spriplasma23","GCF_003339775.1_ASM333977v1_genomic.fna.cmb":"spriplasma24","GCF_001267155.1_ASM126715v1_genomic.fna.cmb":"spriplasma25","GCF_003363775.1_ASM336377v1_genomic.fna.cmb":"spriplasma26","GCF_001274875.1_ASM127487v1_genomic.fna.cmb":"spriplasma27","GCF_000953015.1_Candidatus_Methylopumilus_turicensis_MMS-10A-171_genomic.fna.cmb":"candidatus1","GCF_000981505.1_Candidatus_Methylopumilus_planktonicus_MMS-2-53_genomic.fna.cmb":"candidatus2","GCF_000007025.1_ASM702v1_genomic.fna.cmb":"rickettsia1","GCF_000283995.1_ASM28399v1_genomic.fna.cmb":"rickettsia2","GCF_000008045.1_ASM804v1_genomic.fna.cmb":"rickettsia3","GCF_000284055.1_ASM28405v1_genomic.fna.cmb":"rickettsia4","GCF_000012145.1_ASM1214v1_genomic.fna.cmb":"rickettsia5","GCF_000284075.1_ASM28407v1_genomic.fna.cmb":"rickettsia6","GCF_000012385.1_ASM1238v1_genomic.fna.cmb":"rickettsia7","GCF_000284155.1_ASM28415v1_genomic.fna.cmb":"rickettsia8","GCF_000014345.1_ASM1434v1_genomic.fna.cmb":"rickettsia9","GCF_000284175.1_ASM28417v1_genomic.fna.cmb":"rickettsia10","GCF_000016625.1_ASM1662v1_genomic.fna.cmb":"rickettsia11","GCF_000284195.1_ASM28419v1_genomic.fna.cmb":"rickettsia12","GCF_000017445.4_ASM1744v3_genomic.fna.cmb":"rickettsia13","GCF_000363905.1_ASM36390v1_genomic.fna.cmb":"rickettsia14","GCF_000018205.1_ASM1820v1_genomic.fna.cmb":"rickettsia15","GCF_000367405.1_ASM36740v1_genomic.fna.cmb":"rickettsia16","GCF_000018225.1_ASM1822v1_genomic.fna.cmb":"rickettsia17","GCF_000499665.2_RMONA_1_genomic.fna.cmb":"rickettsia18","GCF_000018245.1_ASM1824v1_genomic.fna.cmb":"rickettsia19","GCF_000831525.1_ASM83152v1_genomic.fna.cmb":"rickettsia20","GCF_000021525.1_ASM2152v1_genomic.fna.cmb":"rickettsia21","GCF_000831545.1_ASM83154v1_genomic.fna.cmb":"rickettsia22","GCF_000022785.1_ASM2278v1_genomic.fna.cmb":"rickettsia23","GCF_001273795.1_ASM127379v1_genomic.fna.cmb":"rickettsia24","GCF_000023005.1_ASM2300v1_genomic.fna.cmb":"rickettsia25","GCF_001442475.1_ASM144247v1_genomic.fna.cmb":"rickettsia26","GCF_000195735.1_ASM19573v1_genomic.fna.cmb":"rickettsia27","GCF_001602215.1_ASM160221v1_genomic.fna.cmb":"rickettsia28","GCF_000221205.1_ASM22120v1_genomic.fna.cmb":"rickettsia29","GCF_001950995.1_ASM195099v1_genomic.fna.cmb":"rickettsia30","GCF_000237845.1_ASM23784v1_genomic.fna.cmb":"rickettsia31","GCF_001951015.1_ASM195101v1_genomic.fna.cmb":"rickettsia32","GCF_000252365.1_ASM25236v1_genomic.fna.cmb":"rickettsia33","GCF_002078315.1_ASM207831v1_genomic.fna.cmb":"rickettsia34","GCF_000277165.1_ASM27716v1_genomic.fna.cmb":"rickettsia35","GCF_002078335.1_ASM207833v1_genomic.fna.cmb":"rickettsia36","GCF_000277185.1_ASM27718v1_genomic.fna.cmb":"rickettsia37","GCF_002285905.1_ASM228590v1_genomic.fna.cmb":"rickettsia38","GCF_000277205.1_ASM27720v1_genomic.fna.cmb":"rickettsia39","GCF_002356695.1_ASM235669v1_genomic.fna.cmb":"rickettsia40","GCF_000277245.1_ASM27724v1_genomic.fna.cmb":"rickettsia41","GCF_002357115.1_ASM235711v1_genomic.fna.cmb":"rickettsia42","GCF_000277265.1_ASM27726v1_genomic.fna.cmb":"rickettsia43","GCF_002357135.1_ASM235713v1_genomic.fna.cmb":"rickettsia44","GCF_000277285.1_ASM27728v1_genomic.fna.cmb":"rickettsia45","GCF_002357155.1_ASM235715v1_genomic.fna.cmb":"rickettsia46","GCF_000277305.1_ASM27730v1_genomic.fna.cmb":"rickettsia47","GCF_002357175.1_ASM235717v1_genomic.fna.cmb":"rickettsia48","GCF_000283595.1_ASM28359v1_genomic.fna.cmb":"rickettsia49","GCF_002357195.1_ASM235719v1_genomic.fna.cmb":"rickettsia50","GCF_000283775.1_ASM28377v1_genomic.fna.cmb":"rickettsia51","GCF_002357215.1_ASM235721v1_genomic.fna.cmb":"rickettsia52","GCF_000283795.1_ASM28379v1_genomic.fna.cmb":"rickettsia53","GCF_002357235.1_ASM235723v1_genomic.fna.cmb":"rickettsia54","GCF_000283815.1_ASM28381v1_genomic.fna.cmb":"rickettsia55","GCF_002357255.1_ASM235725v1_genomic.fna.cmb":"rickettsia56","GCF_000283915.1_ASM28391v1_genomic.fna.cmb":"rickettsia57","GCF_002357275.1_ASM235727v1_genomic.fna.cmb":"rickettsia58","GCF_000283935.1_ASM28393v1_genomic.fna.cmb":"rickettsia59","GCF_002357295.1_ASM235729v1_genomic.fna.cmb":"rickettsia60","GCF_000283955.1_ASM28395v1_genomic.fna.cmb":"rickettsia61","NC_002978.6":"wolbachia1","NC_006833.1":"wolbachia2","NC_018267.1":"wolbachia3","NZ_AP013028.1":"wolbachia4","NZ_DS996944.1":"wolbachia5","NZ_MJMG01000001.1":"wolbachia6"}


class CommandLine():
		"""
		SCAFFOLD CREATED BY DAVID BERNICK
		Handle the command line, usage and help requests.
		CommandLine uses argparse, now standard in 2.7 and beyond.
		it implements a standard command line argument parser with various argument options,
		a standard usage and help, and an error termination mechanism do-usage_and_die.
		attributes:
		all arguments received from the commandline using .add_argument will be
		avalable within the .args attribute of object instantiated from CommandLine.
		For example, if myCommandLine is an object of the class, and requiredbool was
		set as an option using add_argument, then myCommandLine.args.requiredbool will
		name that option.
		"""

		def __init__(self, inOpts=None):
				'''
				CommandLine constructor.
				Implements a parser to interpret the command line argv string using argparse.
				'''

				import argparse
				self.parser = argparse.ArgumentParser(description ='',
				epilog = '',
				add_help = True, #default is True
				prefix_chars = '-',
#       usage = '%(prog)s [options] -option1[default] <input >output'
				)

				self.parser.add_argument('-b', '--blast',required=True,type=str,help='path to blast output file [required]')
				self.parser.add_argument('-f', '--fasta',required=True,type=str,help='query fasta file [required]')
				#self.parser.add_argument('-s', '--spots', required=True, type=int, help="number of total spots [required]")
				self.parser.add_argument('-e','--end',required=True,type=str,help="pair or single end reads [required]")
				if inOpts is None:
						self.args = vars(self.parser.parse_args())
				else:
						self.args = vars(self.parser.parse_args(inOpts))

class Read:
		"""
		A class for reads.
		"""
		def __init__(self,name,subject_hit,length,sstart,send,e_score):
				self.name = name
				self.length = int(length)
				self.subject_hit = subject_hit
				# self.host = host
				self.sstart = sstart
				self.send = send
				self.e_score = e_score



def parseQueryFasta(fasta_file):
	download_num = 0 
	with open(fasta_file,"r") as f: 
			for line in f: 
					line = line.rstrip()
					if line: 
						if line[0] == ">":
								download_num += 1 
	return download_num 

def prune(read_objects):
		"""
		Takes in all the reads and a reference and prunes the best reads that hit the reference.
		Best reads are unique, reads that hit multiple places are ranked by e-value, the lowest e-value wins.
		Returns a list of reads where each object is a best match for that read.
		"""
		import numpy as np
		best_reads = []
		read_name_reads = {}
		for read in read_objects:
				if read.name in read_name_reads:
						read_name_reads[read.name].append(read)
				else:
						read_name_reads[read.name] = [read]
		for read_name, read_objects in read_name_reads.items():
				e_scores = [read.e_score for read in read_objects]
				i = np.argmin(e_scores)
				best_reads.append(read_objects[i])
		return best_reads

def parseBlastOutput(blast_path):
		"""
		This function takes a path to an srr blast output file and returns the best reference and best unique reads.
		The reads are returned as a list of read objects.
		The best reference is determined by the reference with the most unique and significant blast hits.
		"""
		#unpruned_read_objects = {}
		#ref_pruned_reads = {}

		unpruned_read_objects = {key:[] for key in COMMON_NAME.keys()}
		ref_pruned_reads = {key:[] for key in COMMON_NAME.keys()}
		with open(blast_path,"r") as f:
				for line in f:

						line = line.rstrip()
						line = line.rsplit()
						# print(line, file=sys.stderr,flush=True)
						if len(line) > 1:
								read_name = line[0]
								subject_hit = line[1]
								length = int(line[3])
								# sstart = int(line[6])
								# send = int(line[7])
								sstart = int(line[8])
								send = int(line[9])
								e_score = float(line[10])

								# CREATE A READ OBJECT FOR EACH OF THESE SIGNIFICANT HITS TO WOLBACHIA ENDOSYMBIONT.
								# IF A READ HITS THE SAME SUBJECT MORE THAN ONCE,
								# SAVE ONLY THE MOST SIGNIFICANT HIT (LOWEST E-SCORE).
								if e_score < 1e-10 and length > 40:
										# if subject_hit in ENDOSYMBIONT_IDS:
										# wol_host = ENDOSYMBIONT_IDS[subject_hit]
										current_read = Read(read_name,subject_hit,length,sstart,send,e_score)
										if subject_hit in unpruned_read_objects:
												unpruned_read_objects[subject_hit].append(current_read)
										else:
												unpruned_read_objects[subject_hit] = [current_read]
		if len(unpruned_read_objects) > 0:
				for ref in unpruned_read_objects.keys():
						pruned_reads_ref = prune(unpruned_read_objects[ref])
						ref_pruned_reads[ref] = pruned_reads_ref

				return unpruned_read_objects, ref_pruned_reads
		else:
				return None, None


def coverageVariance(reads,reference):
		import numpy as np
		from math import sqrt

		binSizeBp = 5000
		genome_size = WOL_GENOME_SIZES[reference]
		bin_margins = np.arange(0,genome_size,binSizeBp)
		bins =  [0] *  len(bin_margins) # List of empty list to add reads to bins

		for read in reads:
				# WHERE DOES THE SHIFTED READ START
				# for i in range(len(ordering[reference])):
						# if read.subject_hit == ordering[reference][i]:
								# to_add = sum([contig_offset[contig] for contig in ordering[reference][:i]])
								# shifted_start = to_add + read.sstart
								# bin_index = int(shifted_start/binSizeBp)
				bin_index = int(read.sstart/binSizeBp)
				bins[bin_index] += 1

		variance = np.var(bins) # CALCULATE THE VARIANCE OF NUMBER OF READS WITHIN EACH BIN
		normalized_variance = sqrt(variance) / np.mean(bins) # NORMALIZE BY THE AVERAGE NUMBER OF READS ACROSS BINS
		return bins, variance, normalized_variance


def main(myCommandLine=None):
		if myCommandLine is None:
				myCommandLine = CommandLine()
		else:
				myCommandLine = CommandLine()

		# SET VARIABLES FROM COMMANDLINE
		blast_path = myCommandLine.args['blast']
		fasta_file = myCommandLine.args['fasta']
		#spots = myCommandLine.args['spots']
		end = myCommandLine.args['end'].upper()
		target = blast_path.rsplit(".")[0]
		min_read_count = 1

		######################################################
		# Print script name and input variables to stderr.
		print(datetime.datetime.now(), file=sys.stderr,flush=True)
		print("#######################################################################", file=sys.stderr,flush=True)
		print("RUNNING:", sys.argv[0],file=sys.stderr,sep="\t",flush=True)
		for variable, my_variable in myCommandLine.args.items():
				print(variable.upper(), my_variable, sep="\t", file=sys.stderr,flush=True)
		print("#######################################################################", file=sys.stderr,flush=True)
		print(file=sys.stderr,flush=True)
		######################################################

		############### Parse BLAST output ###############
		print("COUNTING NUMBER OF READS SAMPLED...", fasta_file, file=sys.stderr, flush=True)
		download_num = parseQueryFasta(fasta_file)
		if download_num == 0:
			print("ERROR: 0 READS IN QUERY FASTA FILE.", fasta_file, file=sys.stderr,flush=True)
			sys.exit(0)
		else:
			print("NUMBER OF SAMPLED READS:", download_num, file=sys.stderr,flush=True)
		######################################################

		############### Parse BLAST output ###############
		print("PARSING BLAST OUTPUT...", blast_path, file=sys.stderr, flush=True)
		unpruned_reads, ref_pruned_reads = parseBlastOutput(blast_path)
		######################################################
		
		print(file=sys.stderr,flush=True)
		# Print Unpruned read counts
		print("UNPRUNED READ COUNTS:", file=sys.stderr, flush=True)
		if unpruned_reads:
			for ref, read_objects in unpruned_reads.items():
					print(ref, len(read_objects), file=sys.stderr, sep="\t", flush=True)
		else:
			print("WARNING: NO BLAST HITS.", file=sys.stderr, sep="\t", flush=True)
			sys.exit(0)

		print(file=sys.stderr,flush=True)

		# Print Pruned read counts
		print("PRUNED READ COUNTS:", file=sys.stderr, flush=True)
		if ref_pruned_reads:
			for ref, read_objects in ref_pruned_reads.items():
					print(ref, len(read_objects), file=sys.stderr, sep="\t", flush=True)
		else:
			print("ERROR: NO SIGNIFICANT BLAST HITS.", file=sys.stderr, sep="\t", flush=True)
			sys.exit(0)

		print(file=sys.stderr,flush=True)
		print("#######################################################################", file=sys.stderr,flush=True)
		######################################################
		
		print("SUMMARY STATISTIC CALCULATIONS", file=sys.stderr, flush=True)
		# Calculate summary statistics for each reference.
		for ref, read_objects in ref_pruned_reads.items():
				human_readable_ref = COMMON_NAME[ref]
				number_blast_hit_bases = sum([read_object.length for read_object in read_objects])
				size_of_wol_genome = WOL_GENOME_SIZES[ref]
				# Estimate Coverage
				if end == "PAIRED": 
					#end_spots = spots*2
					end_spots = download_num 
				elif end == "SINGLE":
					#end_spots = spots
					end_spots = download_num

				print(ref,human_readable_ref,number_blast_hit_bases,end_spots,download_num,size_of_wol_genome,file=sys.stderr,flush=True,sep="\t")
				estimated_coverage = ( (number_blast_hit_bases * end_spots ) / download_num ) / size_of_wol_genome
				# Estimate Variance
				bins, variance, variance_coefficient = coverageVariance(read_objects,ref)
				percent_bins = len([bin for bin in bins if bin >= min_read_count]) / len(bins)
				print(target, human_readable_ref, ref, end_spots ,download_num, number_blast_hit_bases, size_of_wol_genome,estimated_coverage, variance_coefficient, percent_bins, bins, len(read_objects), sep = "\t", flush=True)
		print("Done.", target, sep="\t", file=sys.stderr,flush=True)

if __name__ == "__main__":
		main()
