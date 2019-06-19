import sys
import os


def countReads(input_reads):
	#read_names = []
	#headers = {}
	#reads = {}
	count = 0 
	with open(input_reads, "r") as f:
		for line in f:
			if line: 	
				line = line.rstrip()
				if line[0] == ">":
					count += 1
					#header = line
					#name = header.rsplit(" ")[0][1:] # name is with the .1 or .2 at the end.		 
					#read_names.append(".".join(name.rsplit(".")[:-1])) # read name is without
					#headers[name] = header 
				#elif line[0] == "@":
				#	header = line
				#	name = header.rsplit(" ")[0][1:] # name is with the .1 or .2 at the end.     
				#	read_names.append(".".join(name.rsplit(".")[:-1])) # read name is without
				#	headers[name] = header
				else:
					pass
					#if name in reads:
					#	reads[name].append(line)
					#else:
					#	reads[name] = [line]
	#return read_names, headers, reads		
	return count 

def printReads(input_reads, count, n): 
	import random 
	probability = n / count 
	
	read_names = []
	headers = {}
	reads = {}
	with open(input_reads,"r") as f: 
		for line in f: 
			if line:
				line = line.rstrip()
				if line[0] == ">":
					print_me = False
					if random.random() < probability:
						print(line, flush=True) 
						print_me = True 
				else: 
					if print_me == True: 
						print(line, flush=True) 

def subsampleReads(read_names, headers, reads, n):
	import random
	subsampled_reads = {}
	subsampled_read_names = random.sample(read_names, n)
	for my_name in subsampled_read_names:
		forward_read = my_name + ".1"
		reverse_read = my_name + ".2"
		print(headers[forward_read],flush=True)
		print("\n".join(reads[forward_read]),flush=True)
		print(headers[reverse_read],flush=True)
		print("\n".join(reads[reverse_read]),flush=True)
	
def main():
	input_reads = sys.argv[1]
	n = 2000000
	count = countReads(input_reads)
	printReads(input_reads,count,n)
	#subsampleReads(read_names, headers, reads, n)

if __name__ == "__main__":
	main()
