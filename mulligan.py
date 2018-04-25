import numpy as np
import random
from numpy.random import uniform
import math
from scipy import optimize
from scipy.special import comb 
from scipy.stats import binom 
from scipy.stats import uniform
import sys
import matplotlib
import cPickle #this is to save stuff
import argparse
import warnings
import scipy
from scipy.sparse.linalg import expm_multiply

#deltas are an array [delta1, delta2, delta3, ..., delta9] which are the parameters of the distribution
#parser = argparse.ArgumentParser()
#parser.add_argument("-f", dest = "filename")
#parser.add_argument("-d", type = float, dest = "distance")
#args = parser.parse_args()
#filename = args.filename
#distance = args.distance

def find_p(n, x, t1, t2):
	h = []
	for i in list(range(n+1)):
		h.append((x**(i))*((1-x)**(n-i))) 
	h = np.asarray(h)
	Q_t2 = np.zeros((n+1, n+1))
	for i in range(n+1): 
		Q_t2[i][i] = -i*(n-i)
		if i >= 1: Q_t2[i][i-1] = 0.5*i*(i-1)
		if i < n: Q_t2[i][i+1] = 0.5*(n-i)*(n-i-1)
	Q_t2 = Q_t2 * t2 

	Q_t1 = np.zeros((n+1, n+1))
	for i in range(n+1): 
		Q_t1[i][i] = -i*(n-i+1)
		if i < n: Q_t1[i][i+1] = 0.5*(n-i+1)*(n-i)
		if i >= 1: Q_t1[i][i-1] = 0.5*i*(i-1)
	Q_t1 = Q_t1 * t1

	Q_t1h = scipy.sparse.linalg.expm_multiply(Q_t1, h)
	p = scipy.sparse.linalg.expm_multiply(Q_t2, Q_t1h) 
		
	return p 

def populate_prob_table(allele_frequencies, t1, t2):
	table = np.zeros((len(allele_frequencies),9,9))
	for i in range(len(allele_frequencies)): 
		x = allele_frequencies[i]
		n1 = find_p(1, x, t1, t2)
		n2 = find_p(2, x, t1, t2)
		n3 = find_p(3, x, t1, t2)
		n4 = find_p(4, x, t1, t2)

		table[i,0] = (n1[0], n2[0], n2[0], n3[0], n2[0], n3[0], n2[0], n3[0], n4[0])
		table[i,1] = (0, n2[1], 0, n3[2], 0, n3[1], 0, 0, n4[2])
		table[i,2] = (0, 0, n2[1], 2*n3[1], 0, 0, 0, n3[1], 2*n4[1])
		table[i,3] = (0, 0, 0, 0, n2[1], 2*n3[1], 0, n3[1], 2*n4[1])
		table[i,4] = (0, 0, 0, 0, 0, 0, 2*n2[1], n3[2] + n3[1], 4*n4[2])
		table[i,5] = (0, 0, 0, 0, n2[1], 2*n3[1], 0, n3[2], 2*n4[3])
		table[i,6] = (0, 0, n2[1], 2*n3[2], 0, 0, 0, n3[2], 2*n4[3])
		table[i,7] = (0, n2[1], 0, n3[1], 0, n3[2], 0, 0, n4[2])
		table[i,8] = (n1[1], n2[2], n2[2], n3[3], n2[2], n3[3], n2[2], n3[3], n4[4])
	#rows are IBS, columns are IBD
#		print table[i]
#		raw_input()
	return table
	
def create_data_given_IBS(N_1, K_1, N_2, K_2, err1, err2, allele_frequencies):
	data_given_IBS = []
	for i in range(len(allele_frequencies)):
		num_of_sites = len(N_1[i]) #( == len(K_1[i]) == lenN_2[i]) = len(K_2[i]) 
	#	print "num of sites =", num_of_sites	
		data_given_IBS.append(np.zeros((num_of_sites,9)))
	#	print "N1 length =", len(N_1[i]), "K1 length =", len(K_1[i])
		pmf_K1_N1_err = binom.pmf(K_1[i], N_1[i], err1)
		pmf_K2_N2_err = binom.pmf(K_2[i], N_2[i], err2)
		pmf_K1_N1_05 = binom.pmf(K_1[i], N_1[i], 0.5)
		pmf_K2_N2_05 = binom.pmf(K_2[i], N_2[i], 0.5)
		pmf_K1_N1_minus_err = binom.pmf(K_1[i], N_1[i], 1-err1)
		pmf_K2_N2_minus_err = binom.pmf(K_2[i], N_2[i], 1-err2)
		data_given_IBS[i] = np.transpose((pmf_K1_N1_err * pmf_K2_N2_err,
				pmf_K1_N1_err * pmf_K2_N2_minus_err,
				pmf_K1_N1_err * pmf_K2_N2_05,
				pmf_K1_N1_05 * pmf_K2_N2_err,
				pmf_K1_N1_05 * pmf_K2_N2_05,
				pmf_K1_N1_05 * pmf_K2_N2_minus_err,
				pmf_K1_N1_minus_err * pmf_K2_N2_05,
				pmf_K1_N1_minus_err * pmf_K2_N2_err,
				pmf_K1_N1_minus_err * pmf_K2_N2_minus_err)) 
#		print "data_given_IBS[i]", data_given_IBS[i][0]
#		raw_input()
	return data_given_IBS 
	
def inbreeding_likelihood(N_1, K_1, N_2, K_2, err1, err2, t1, t2, allele_frequencies, delta):
#	should be the modern allele frequency
	delta_sum = sum(delta)
	if delta_sum > 1: 
		likelihood = -1e300
		return likelihood
	delta = np.append(delta,(1-delta_sum))
	IBS_given_IBD = populate_prob_table(allele_frequencies, t1, t2)
	data_given_IBS = create_data_given_IBS(N_1, K_1, N_2, K_2, err1, err2, allele_frequencies)
#	print data_given_IBS[0][0]
	#calculating prob_IBS and prob_data should happen in the loop over allele frequcies
#	cur_likelihood = 0.
	loglikelihood = 0.
	cur_likelihood = []
	for i in range(len(allele_frequencies)):
		if allele_frequencies[i] == 0 or allele_frequencies[i] == 1: continue
#		with warnings.catch_warnings(record = True):
#			warnings.filterwarnings('error')
#			warnings.simplefilter('default')
#			try:
		prob_IBS = np.dot(IBS_given_IBD[i], delta)
#				print "len prob_IBS=", len(prob_IBS), "len data_given_IBS[i] =", len(data_given_IBS[i])
		prob_data = np.dot(data_given_IBS[i], prob_IBS)
#				print data_given_IBS[i], " ---\n", prob_IBS
#				raw_input()
		cur_likelihood = np.log(prob_data)
#				warnings.warn(Warning())
#		print cur_likelihood
#			except Warning:
		which_zero = np.where(prob_data == 0)[0]
		if len(which_zero) >= 1:
			which_zero = which_zero[:10]
#			print "PROB IBS GIVEN IBD"
#			print IBS_given_IBD[i]
#			print "PROB IBS"
#			print prob_IBS 
#			print "NUMBER OF ZEROS AND NUMBER OF DATA"
#			print len(which_zero), len(prob_data)
#			print "DATA"
#			print which_zero, N_1[i][which_zero], K_1[i][which_zero], N_2[i][which_zero], K_2[i][which_zero]
#			print "PROB DATA GIVEN IBS"
#			print data_given_IBS[i][which_zero]
#			print "PROB DATA"
#			print prob_data[which_zero]
#			print allele_frequencies[i]
#			raw_input()
#				print "Dot product of 0", delta 
#				return -np.inf
#			finally:
		loglikelihood += np.sum(cur_likelihood)
#	print delta, t1, t2, err1, err2, loglikelihood
	return loglikelihood

#TO SAVE INTO A PICKLE
#pickle_file = open("LnL.pickle","w")
#cPickle.dump(LnL,pickle_file)
#pickle_file.close()

#TO LOAD
#pickle_file = open("LnL.pickle")
#LnL = cPickle.load(pickle_file)
#pickle_file.close()



def main(delta7, delta8, N_1, K_1, N_2, K_2, allele_frequencies):

	t = [uniform.rvs(loc = 1e-7, scale = 2, size = 2)]
	e = [uniform.rvs(loc = 1e-6, scale = 0.5, size = 2)]
	x0 = np.append(t, e)

#pars = [tl, t2, err1, err2]

#	opt = optimize.fmin_l_bfgs_b(func = lambda pars: -inbreeding_likelihood(N_1, K_1, N_2, K_2, pars[2], pars[3], pars[0], pars[1], allele_frequencies, [0,0,0,0,0,0,delta7,delta8]), x0 = x0, approx_grad = True, epsilon = 1e-5, bounds =  [[1e-6,2]]*2 + [[1e-6,.5]]*2, factr = 1, pgtol = 1e-12)	 
	opt = optimize.fmin_l_bfgs_b(func = lambda pars: -inbreeding_likelihood(N_1, K_1, N_2, K_2, pars[0], pars[1], 0, 0, allele_frequencies, [0,0,0,0,0,0,delta7,delta8]), x0 = x0, approx_grad = True, epsilon = 1e-5, bounds =  [[1e-6,2]]*2 + [[1e-6,.5]]*2, factr = 1, pgtol = 1e-12)	 
#opt = optimize.fmin_tnc(func = lambda pars: -inbreeding_likelihood(N_1, K_1, N_2,K_2, pars[4], pars[5], pars[2], pars[3], allele_frequencies, [0, 0, 0, 0, 0, 0, pars[0], pars[1]]), x0 = x0, approx_grad = True, bounds = ([[1e-6,1]]*2 + [[1e-6,2]]*2 + [[1e-6,0.5]]*2),  epsilon = 1e-8, maxfun = 1000)
#	opt = optimize.fmin_cobyla(func = lambda pars: -inbreeding_likelihood(N_1, K_1, N_2, K_2, pars[4], pars[5], pars[2], pars[3], allele_frequencies, [0,0,0,0,0,0,pars[0],pars[1]]), x0 = x0, cons = [constr1, constr2, constr3, constr4, constr5, constr6, constr7, constr8, constr9], rhobeg = 0.05, maxfun=10000) 
 
	print opt
	return opt
