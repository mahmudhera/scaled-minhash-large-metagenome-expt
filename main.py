import screed
import subprocess
import mmh3
import numpy as np
import pandas as pd

def get_kmers_in_file(filename, k):
	# run jellyfish
	cmd = "jellyfish count -m " + str(k) + " -s 2G -t 32 -o tmp -C -L 1 -U 99999999 " + filename
	args = cmd.split(' ')
	subprocess.call(args)
	cmd = "jellyfish dump -L 1 -U 99999999 -o tmp-dump -c tmp"
	args = cmd.split(' ')
	subprocess.call(args)
	# open output file
	df = pd.read_csv('tmp-dump', delimiter=' ', header=None)
	list_kmers = df.iloc[:,0].tolist()
	return list_kmers
	# construct set
	# return the set
	

# c is float, 0 < c < 1
def generate_c_percent_of_file(c, genome_filename, out_filename):
    subprocess.call(['rm', out_filename])
    length = 0
    with screed.open(genome_filename) as f:
        for record in f:
            length += len(record.sequence)
    required_length = int(length * c)
    with screed.open(genome_filename) as f:
        for record in f:
            if len(record.sequence) > required_length:
                small_str = record.sequence[:required_length]
                break
    f2 = open(out_filename, 'w')
    f2.write('>small_seq\n')
    f2.write(small_str)
    f2.write('\n')
    f2.close()
	
	
if __name__ == "__main__":
	mg_filename = "reads_file.fastq"
	g_filename = 'staphylococcus.fasta'
	smallg_filename = 'temp.fasta'
	k = 21
	containment_ranges = [0.01] + [0.1*i for i in range(1, 10)] + [0.99]
	
	a1 = get_kmers_in_file(mg_filename, k)
	b = get_kmers_in_file(g_filename, k)
	
	print(len(a1), a[0])
	print(len(b), b[0])
	
	print("---")
	
	for C in containment_ranges:
		extract_part_of_genome(C, g_filename, smallg_filename)
		a2 = get_kmers_in_file(smallg_filename, k)
		print(len(a2), a2[0])