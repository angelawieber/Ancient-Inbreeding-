import numpy as np
import math

#deltas are an array [delta1, delta2, delta3, ..., delta9] which are the parameters of the distribution
#genotype data should be a list of numbers 0 through 4 corresponding to the appropriate genotype configurations. Each entry is a SINGLE SITE
#allele_frequencies is an array of the frequencies of one of the alleles at every site
genotype_data = [0, 1, 2, 3, 4]
def inbreeding_likelihood(delta, genotype_data, allele_frequencies):
	if len(genotype_data) != len(allele_frequencies):
		print "ERROR: genotype data and allele frequencies not the same length"
		return None
	likelihood = 0
	for i in range(len(genotype_data)): 
		cur_likelihood = np.log(inbreeding_likelihood_one_site(delta, genotype_data[i], allele_frequencies[i]))
		print cur_likelihood
		likelihood += cur_likelihood
	return likelihood

def populate_prob_table(allele_freq):
	table = np.zeros((9,9))
	p = allele_freq
	print p
	table[0] = (1-p, (1-p)**2, (1-p)**2, (1-p)**3, (1-p)**2, (1-p)**3, (1-p)**2, (1-p)**3, (1-p)**4)
	table[1] = (0, p*(1-p), 0, (1-p)*(p**2), 0, p*(1-p)**2, 0, 0, (p**2)*(1-p)**2)
	table[2] = (0, 0, p*(1-p), 2*p*(1-p)**2, 0, 0, 0, p*(1-p)**2, 2*p*(1-p)**3)
	table[3] = (0, 0, 0, 0, p*(1-p), 2*p*(1-p)**2, 0, p*(1-p)**2, 2*p*(1-p)**3)
 	table[4] = (0, 0, 0, 0, 0, 0, 2*p*(1-p), p*(1-p), 4*(p**2)*(1-p)**2)
	table[5] = (0, 0, 0, 0, p*(1-p), 2*(p**2)*(1-p), 0, (p**2)*(1-p), 2*(p**3)*(1-p))
	table[6] = (0, 0, p*(1-p), 2*(p**2)*(1-p), 0, 0, 0, (p**2)*(1-p), 2*(p**3)*(1-p))
	table[7] = (0, p*(1-p), 0, p*(1-p)**2, 0, (p**2)*(1-p), 0, 0, (p**2)*(1-p)**2)
	table[8] = (p, p**2, p**2, p**3, p**2, p**3, p**2, p**3, p**4)

	return table


#delta is still all the deltas
#genotype is just ONE genotype, for one site
#allele_freq is just ONE allel frequency 
def inbreeding_likelihood_one_site(delta, genotype, allele_freq):
	prob_table = populate_prob_table(allele_freq)
	likelihood = np.dot(prob_table[genotype], delta)
	return likelihood

test_like = inbreeding_likelihood([1/9.]*9,[3,1,3,2,0,4],[.5,.1,.3,.7,.8,.2])
print test_like



	
