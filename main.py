import screed
import subprocess
import mmh3
import numpy as np
import pandas as pd


class ScaledMinHash:
	def __init__(self, scale_factor, max_hash_value, initial_set=None):
		if initial_set is None:
			self.hash_set = set()
		else:
			self.hash_set = set(initial_set)
		self.H = max_hash_value
		self.scale_factor = scale_factor
		self.raw_elements = set()
	
	def add_value(self, hash_value):
		if hash_value <= self.H * self.scale_factor:
			self.hash_set.add(hash_value)
		self.raw_elements.add(hash_value)
	
	def add_values(self, hash_values):
		for hash_value in hash_values:
			self.add_value(hash_value)
    
	def remove(self, hash_value):
		self.hash_set -= hash_value
	
	def get_containment(self, smh):
		return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / len(self.hash_set)
		
	def get_scaled_containment(self, smh):
		bf = 1 - (1 - self.scale_factor) ** len(self.raw_elements)
		return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / ( len(self.hash_set) * bf )

	def get_sketch_size(self):
		return len( self.hash_set )
		
def get_hash_from_kmer(kmer, seed=0):
	hash_value = mmh3.hash64(kmer, seed=seed)[0]
	if hash_value < 0:
		hash_value += 2**64
	return hash_value
	
	
def create_scaled_minhash(kmers, seed, scale_facor):
	H = 2**64
	smh1 = ScaledMinHash(scale_facor, H)
	for kmer in kmers:
		h = get_hash_from_kmer(kmer, seed)
		smh1.add_value(h)
	return smh1
	
def add_kmers_in_scaled_minhash(kmers, smh, seed):
	H = 2**64
	for kmer in kmers:
		h = get_hash_from_kmer(kmer, seed)
		smh.add_value(h)
	return smh


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
	#mg_filename = "reads_file.fastq"
	mg_filename = 'SRR492065.contigs.fa'
	g_filename = 'staphylococcus.fasta'
	#g_filename = 'SRR492065.contigs.fa'
	smallg_filename = 'temp.fasta'
	k = 21
	#containment_ranges = [0.01] + [0.1*i for i in range(1, 10)] + [0.99]
	containment_ranges = [0.2]
	scale_facor = 0.0001
	num_runs = 2
	seeds = [i for i in range(num_runs)]
	seed = 1
	
	kmers_in_metagenome = get_kmers_in_file(mg_filename, k)
	kmers_in_genome = get_kmers_in_file(g_filename, k)
	
	sketches_metagenome = {}
	print('generating sketches for metagenome')
	for seed in seeds:
		sketches_metagenome[seed] = create_scaled_minhash(kmers_in_metagenome, seed, scale_facor)
	print('sketch_sizes:')
	print([sketches_metagenome[seed].get_sketch_size() for seed in seeds])
	
	sketch_genome = sketches_metagenome[0]
	
	print('generating sketch for kmers in genome')
	sketch_genome = create_scaled_minhash(kmers_in_genome, seed, scale_facor)
	print('sketch size:')
	print(sketch_genome.get_sketch_size())
	
	for C in containment_ranges:
		generate_c_percent_of_file(C, g_filename, smallg_filename)
		kmers_in_small_portion = get_kmers_in_file(smallg_filename, k)
		
		#write fn to create new smh, then add these new kmers
		sketch_new = ScaledMinHash( sketch_metagenome.scale_factor, sketch_metagenome.H, sketch_metagenome.hash_set )
		add_kmers_in_scaled_minhash(kmers_in_small_portion, sketch_new, 1)
		print('added kmers in small genome in sketch. new sletch size:')
		print(sketch_new.get_sketch_size())
		
		# list = []
		
		# for seed in all_seeds:
			# determine num_union
			# know how many kmers in one of them already
			# get mash jaccard
			# get mash containment
			# create two msh sketches, and get scaled_containment
			# add values in the list
		
		# after loop, output the value and variance with C
		
		unique_kmers_union = set(kmers_in_genome + kmers_in_metagenome + kmers_in_small_portion)
		print('Unique kmers in union:')
		print(len(unique_kmers_union))
		
		print("Scaled containment: C(genome in metagenome):")
		print(sketch_metagenome.get_scaled_containment(sketch_genome))
		
		print("Scaled containment: C(genome in metagenome+small):")
		print(sketch_new.get_scaled_containment(sketch_genome))