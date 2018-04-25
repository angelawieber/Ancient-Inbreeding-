from pysam import VariantFile
import sys
import argparse

parser = argparse.ArgumentParser(description = 'Find IBD modes for 2 individuals.')
parser.add_argument("-v1", action = "store", dest = "vcf1")
parser.add_argument("-v2", dest = "vcf2")
parser.add_argument("-v3", dest = "vcf3")
parser.add_argument("-v4", dest = "vcf4")
parser.add_argument("-v5", dest = "vcf5")
parser.add_argument("-v6", dest = "vcf6")
parser.add_argument("-v7", dest = "vcf7")
parser.add_argument("-v8", dest = "vcf8")
parser.add_argument("-v9", dest = "vcf9")
parser.add_argument("-v10", dest = "vcf10")
parser.add_argument("-v11", dest = "vcf11")
parser.add_argument("-v12", dest = "vcf12")
parser.add_argument("-v13", dest = "vcf13")
parser.add_argument("-v14", dest = "vcf14")
parser.add_argument("-v15", dest = "vcf15")
parser.add_argument("-v16", dest = "vcf16")
parser.add_argument("-v17", dest = "vcf17")
parser.add_argument("-v18", dest = "vcf18")
parser.add_argument("-v19", dest = "vcf19")
parser.add_argument("-v20", dest = "vcf20")
parser.add_argument("-v21", dest = "vcf21")
parser.add_argument("-v22", dest = "vcf22")
parser.add_argument('-o', action = "store", dest = "output_filename")
parser.add_argument('-p1', action = "store", dest = "pop")
#parser.add_argument('-p2', action = "store", dest = "pop2")
parser.add_argument('-i', action = "store", dest = "ind1")
parser.add_argument('-j', action = "store", dest = "ind2")
args = parser.parse_args()
vcf_filename = [args.vcf1, args.vcf2, args.vcf3, args.vcf4, args.vcf5, args.vcf6, args.vcf7, args.vcf8, args.vcf9, args.vcf10, args.vcf11, args.vcf12, args.vcf13, args.vcf14, args.vcf15, args.vcf16, args.vcf17, args.vcf18, args.vcf19, args.vcf20, args.vcf21, args.vcf22]
output_filename = args.output_filename
pop = args.pop
#pop2 = args.pop2
ind1 = args.ind1
ind2 = args.ind2
#print args


f = open("/mnt/sda/data/1KG/samples_and_pops.txt")
counter = 0
outfile = open(output_filename,"a")	

IDS = []
for line in f:
	sub_pop = line.split()[1] 
	if sub_pop == pop: IDS.append(line.split()[0])
print IDS
for vcf in vcf_filename: 
	vcf_in = VariantFile(vcf)
	for rec in vcf_in.fetch():
	#only print biallelic sites (i.e. ones with only 2 alleles)
	#print the DERIVED allele frequency, i.e. the frequency of the non-ancestral allele
	#skp ones where there's no ancestral allele called
		counter += 1 
	#if counter > 2000: break
		if "AA" not in rec.info: continue
		if len(rec.alts) != 1 or len(rec.ref) != 1: continue
	#print rec.chrom, rec.ref, rec.alts[0], rec.pos, rec.info["EUR_AF"][0]
		AAfield = rec.info["AA"]
		AAsplit = AAfield.split("|")
		AA = AAsplit[0].upper() 
		if (AA != 'A' and AA != 'C' and AA != 'T' and AA != 'G'): continue
		alt_freqs = []
		num_alt_alleles = 0.0 #keeps track of how many alternative alleles you've seen so far for this site
		num_good_alleles = 0.0 #keeps track of how man non-missing alleles you've seen so far
		num_none = 0
#		alt_freq = float(rec.info["EUR_AF"][0])
		for ID in IDS:
			genotype = rec.samples[ID].allele_indices
			num_alt_alleles += sum(genotype)
			for geno in genotype:
				 if geno == None: num_none + 1
			num_good_alleles += len(genotype) - num_none
	#count the alleles that aren't None
		freq = num_alt_alleles/num_good_alleles
	#from here it's like you just got he EUR_AF
		if rec.ref != AA: freq = 1 - freq
		alt_freqs.append(freq)
	#print alt_freqs 
		alt_freq = sum(alt_freqs)/len(alt_freqs)
		if alt_freq == 0 or alt_freq == 1: continue
	#get what the derived allele is
		if AA == rec.ref: der = rec.alts[0]
		if len(der) != 1: continue 
		elif AA == rec.alts[0]: der = rec.ref
		switch = {(0,0):(1,1), (0,1):(1,0), (1,0):(0,1), (1,1):(0,0)} 
		GT1 = rec.samples[ind1].allele_indices
		if AA == rec.alts[0]: GT1 = switch[GT1]
		GT2 = rec.samples[ind2].allele_indices
		if AA == rec.alts[0]: GT2 = switch[GT2]
		config = GT1 + GT2
		config_dict = {(0,0,0,0):0, (0,0,1,1):1, (0,0,0,1):2, (0,0,1,0):2, (0,1,0,0):3, (1,0,0,0):3, (0,1,0,1):4, (1,0,0,1):4, (1,0,1,0):4, (0,1,1,0):4,
	 	(1,0,1,1):5, (0,1,1,1):5, (1,1,0,1):6, (1,1,1,0):6, (1,1,0,0):7, (1,1,1,1):8} 	
		config = config_dict[config]
	#map this to the genotype configurations numbers
	#add that as another column to the output.
		outfile.write("%s\t%d\t%s\t%s\t%f\t%d\n"%(rec.chrom,rec.pos,AA,der,alt_freq,config))

