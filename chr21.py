from pysam import VariantFile
vcf_in = VariantFile("/mnt/sda/data/1KG/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
counter = 0
outfile = open("EUR_AF.txt","w")
for rec in vcf_in.fetch():
	#only print biallelic sites (i.e. ones with only 2 alleles)
	#print the DERIVED allele frequency, i.e. the frequency of the non-ancestral allele
	#skp ones where there's no ancestral allele called
	counter += 1 
	if counter > 2000: break
	if "AA" not in rec.info: continue
	if len(rec.alts) != 1 or len(rec.ref) != 1: continue
	#print rec.chrom, rec.ref, rec.alts[0], rec.pos, rec.info["EUR_AF"][0]
	AAfield = rec.info["AA"]
	AAsplit = AAfield.split("|")
	AA = AAsplit[0].upper() #TODO: figure out what lower case means
	if AA == "." or AA == "": continue
	alt_freq = rec.info["EUR_AF"][0] #+ rec.info["AFR_AF"][0] + rec.info["AMR_AF"][0] + rec.info["EAS_AF"][0] + rec.info["SAS_AF"][0] 
	if rec.ref != AA: alt_freq = 1 - alt_freq
	if alt_freq == 0 or alt_freq ==1: continue
	#get what the derived allele is
	if AA == rec.ref: der = rec.alts[0]
	elif AA == rec.alts[0]: der = rec.ref
	switch = {(0,0):(1,1), (0,1):(1,0), (1,0):(0,1), (1,1):(0,0)} 
	GT1 = rec.samples["HG00096"].allele_indices
	if AA == rec.alts[0]: GT1 = switch[GT1]
	GT2 = rec.samples["HG00097"].allele_indices
	if AA == rec.alts[0]: GT2 = switch[GT2]
	config = GT1 + GT2
	config_dict = {(0,0,0,0):0, (0,0,1,1):1, (0,0,0,1):2, (0,0,1,0):2, (0,1,0,0):3, (1,0,0,0):3, (0,1,0,1):4, (1,0,0,1):4, (1,0,1,0):4, (0,1,1,0):4,
 	(1,0,1,1):5, (0,1,1,1):5, (1,1,0,1):6, (1,1,1,0):6, (1,1,0,0):7, (1,1,1,1):8} 	
	config = config_dict[config]
	#map this to the genotype configurations numbers
	#add that as another column to the output.
	outfile.write("%s\t%d\t%s\t%s\t%f\t%s\n"%(rec.chrom,rec.pos,AA,der,alt_freq,config))

