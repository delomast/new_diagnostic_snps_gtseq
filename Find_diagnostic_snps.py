#!/usr/bin/python3
##second step in snp discovery from GT-seq reads pipeline - SPECIFIC TO FINDING HYBRID DIAGNOSTIC MARKERS
##Looks through haplotypes and finds potential diagnostic markers 
#all three input parameters are optional
#default -min_prop is 0.8
#default -min_depth is 10
#default -min_error is 0
# -s1 and -s2 are txt files with sample names in species1 and species2, one name per line
#In terminal type python Find_snps.py -min_prop min_proportion_of_samples_with_haplotypes -s1 file_with_species1_samples 
# 							-s2 file_with_species2_samples -min_depth min_depth_for_calling_haplotypes -min_error proportion of alleles that can be wrong in one group or the other


import sys

def Main():

	samples = {}
	#read in sample names, marker names, and haplotypes
	#all_haplotypes.txt is tab-delimited indiv_name, marker_name, hap_1, depth_1, hap_2, depth_2
	for line in open('all_haplotypes.txt', 'r'):
		line = line.rstrip()
		sep = line.split('\t')
		if sep[0] in samples:
			samples[sep[0]][sep[1]] = sep[2:6]
		else:
			samples[sep[0]] = {sep[1] : sep[2:6]}

			
	markers = []		#get list of markers from one individual (all individuals will have entries for all markers, whether they were typed for them or not)
	for i in samples[list(samples.keys())[0]]:
		markers.append(i)
	
	min_num_haplo = 0.8		#use default of 0.8 if no value is given
	min_maf = 0.05		#use default of 0.05 if no value is given
	min_depth = 10		#use default of 10 if no value is given
	min_error = 0		#use default of 0 if no value is given
	s1_file = 'none'
	s2_file = 'none'
	
	for flag in range(1, len(sys.argv), 1):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == '-min_depth':
			min_depth = int(sys.argv[flag + 1])
		if sys.argv[flag] == '-min_prop':
			min_num_haplo = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-min_error':
			min_error = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-s1':
			s1_file = sys.argv[flag + 1]		
		if sys.argv[flag] == '-s2':
			s2_file = sys.argv[flag + 1]

	if s1_file == 'none':
		print('No species1 sample list given. Exiting.')
		return	
	if s2_file == 'none':
		print('No species2 sample list given. Exiting.')
		return
				
	
	#read in s1 and s2 sample names
	s1_names, s2_names = [], []
	with open(s1_file, 'r') as file_in:
		for line in file_in:
			s1_names.append(line.rstrip())
	with open(s2_file, 'r') as file_in:
		for line in file_in:
			s2_names.append(line.rstrip())		
			
	print('Loaded ', str(len(s1_names)), ' samples in species 1 and ', str(len(s2_names)), ' samples in species 2\n\n')

	print('Finding diagnostic markers with reads in at least ', min_num_haplo, ' proportion of samples\n\t a minimum error rate of ', min_error, '\n\t and using a minimum depth of ', min_depth)

	min_haplo_1 = round(min_num_haplo*len(s1_names)*2,0)
	min_haplo_2 = round(min_num_haplo*len(s2_names)*2,0)
	
	poten_out = open('potential_diag_snps.txt', 'w')
	poten_out.write('locus\tcol\terror_sp1\terror_sp2\n')
	for marker in markers:	#for all markers
		#get all haplotypes for the marker
		all_haplotypes_sp1, all_haplotypes_sp2 = [], []
		for indiv in s1_names:	#and for all samples in s1
			if samples[indiv][marker][0][0] != 'N' and int(samples[indiv][marker][1]) >= min_depth: #check that there is an actual read and that depth is above min_depth
				all_haplotypes_sp1.append(samples[indiv][marker][0])		
			if samples[indiv][marker][2][0] != 'N' and int(samples[indiv][marker][3]) >= min_depth: 
				all_haplotypes_sp1.append(samples[indiv][marker][2])

		for indiv in s2_names:	#and for all samples in s2
			if samples[indiv][marker][0][0] != 'N' and int(samples[indiv][marker][1]) >= min_depth: #check that there is an actual read and that depth is above min_depth
				all_haplotypes_sp2.append(samples[indiv][marker][0])		
			if samples[indiv][marker][2][0] != 'N' and int(samples[indiv][marker][3]) >= min_depth: 
				all_haplotypes_sp2.append(samples[indiv][marker][2])

		#discard markers that didn't have haplotypes in minimum number of individuals
		if len(all_haplotypes_sp1) < min_haplo_1 or len(all_haplotypes_sp2) < min_haplo_2:
			continue
		#find snps 
		max_length = 0
		for i in all_haplotypes_sp1 + all_haplotypes_sp2:
			if len(i) > max_length:
				max_length = len(i)		#find maximum length of haplotypes
		for i in range(0, max_length, 1):
			unique_sp1 = {}					#build dictionary of alleles and their frequency at a given position in the haplotype
			unique_sp2 = {}					#build dictionary of alleles and their frequency at a given position in the haplotype
			for haplo in all_haplotypes_sp1:
				if len(haplo) < (i + 1):		#make sure current haplotype is long enough to have a base
					continue
				if haplo[i] == 'N':		#prevent no call from being counted as an allele
					continue
				if haplo[i] in unique_sp1:
					unique_sp1[haplo[i]] += 1
				else:
					unique_sp1[haplo[i]] = 1
			for haplo in all_haplotypes_sp2:
				if len(haplo) < (i + 1):		#make sure current haplotype is long enough to have a base
					continue
				if haplo[i] == 'N':		#prevent no call from being counted as an allele
					continue
				if haplo[i] in unique_sp2:
					unique_sp2[haplo[i]] += 1
				else:
					unique_sp2[haplo[i]] = 1
			
			#make sure both species have reads at this position
			if len(unique_sp1) == 0 or len(unique_sp2) == 0:
				continue
			#determine if diagnostic
			### express allele counts as proportions
			total_sp1, total_sp2 = 0, 0
			for j in unique_sp1:
				total_sp1 += unique_sp1[j]
			for j in unique_sp2:
				total_sp2 += unique_sp2[j]
			for j in unique_sp1:
				unique_sp1[j] = unique_sp1[j] / total_sp1
			for j in unique_sp2:
				unique_sp2[j] = unique_sp2[j] / total_sp2
			### determine overlap
			overlap_sp1, overlap_sp2 = 0.0, 0.0
			for j in unique_sp1:
				if j in unique_sp2:
					overlap_sp1 += unique_sp1[j]
			for j in unique_sp2:
				if j in unique_sp1:
					overlap_sp2 += unique_sp2[j]
			#output snp and error if it passes error cutoff
			if overlap_sp1 <= min_error and overlap_sp2 <= min_error:
				poten_out.write(marker + '\t' + str(i) + '\t' + str(overlap_sp1) + '\t' + str(overlap_sp2) + '\n')
	
	poten_out.close()

Main()


