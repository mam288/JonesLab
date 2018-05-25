##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:28:50 2017

@author: meisha
"""
import argparse
from pyfasta import Fasta
import time

import sys
sys.path.insert(0,'/nas/longleaf/home/mam288/.local/lib/python3.5/site-packages') # path to python 3


import pandas as pd
import datetime
import os, sys
from itertools import chain, combinations
start_time = time.time()
now = datetime.datetime.now()

def safe_div(x,y):
#  Allows for division by 0. Returns x/y if y != 0 and 0 if y == 0 

    if y==0: return 0
    return x/y

def all_subsets(input_list):
#  Returns a list of all possible combinations of the alleles entered as a list.
#  
#  Input: 
#      input_list (list of alleles)
#  Output: 
#      All possible combinations of the alleles in input_list (list)

    return chain(*map(lambda x: combinations(input_list, x), range(0, len(input_list)+1)))

def all_empty_lists(variant_list):
#    Returns True if all of the sublists in the list entered are empty, otherwise returns False
#    
#    Input: 
#        variant_list - list of variants
#    Output: 
#        boolean - True if all variants in list_var are empty lists, else False

    for l in variant_list:
        if l != []:
            return False
    return True

def is_equivalent(parent_variants, offspring_variants):
#    Returns True if any of the alleles created using the offpsring_variants match any of the alleles created using the parent_variants
#    Input:
#        parent_variants (tuple) - parent variants for comparison
#        offspring_variants (tuple) - offspring variants for comparison
#    Output:
#        boolean - Flase if any of the offpsring alleles do not match any of the parent alleles

    c = 0
    for subset in all_subsets(offspring_variants): # get all subsets of the allele list so all combinations can be compared against parent alleles
        if subset != ():
            if c == 0:
                offspring_variants = (offspring_variants,) + (subset,)
            else:
                offspring_variants = offspring_variants + subset
            c += 1
        
    #If all parent and offspring alleles are empty return True and return False if the parent allele list
    #is empty but the offspring is not. If the offspring is empty but parents are not return 'Empty Allele List' '''
    if all_empty_lists(parent_variants) and all_empty_lists(offspring_variants):
        return True
    if (all_empty_lists(parent_variants) and not all_empty_lists(offspring_variants)):
        return False
    if (all_empty_lists(offspring_variants) and not all_empty_lists(parent_variants)):
        return 'Empty Allele List'
    
    for variant_set in offspring_variants: # go through each set of alleles in the offspring_variants

        match_list = []
        
        for par_var in parent_variants:
            match_list += [is_equivalent_compare_variants(par_var,variant_set)]
        for var in variant_set:
            match_list += [is_equivalent_compare_variants(var,parent_variants)]
        return True in match_list # return True if any of the allleles match
            
def is_equivalent_compare_variants(variant,variant_list_for_comparison):
#    Returns True if any of the variants in variant_list_for_comparison match the variant. Variants: [variant start position, 
#    reference allele, alternate allele, microsatellite start position, microsatellite sequence, average coverage]. 
#    Input: 
#        variant (list) - variant 
#        variant_list_for_comparison (list of variants) - list of variants to compare to
#    Output:
#        Boolean value indicating whether any of the variants in variant_list_for_comparison match the inputted variant.

    if variant_list_for_comparison == ():
        if variant == ():
            return True # if both variant and variant_list_for_comparison are empty return True
        else:
            return False # if a2 is empty but variant is not

    new_seq_list = [] # list of reconstructed alleles

    c = 0 

    for var in (variant,) + variant_list_for_comparison:

        # minimum difference in length between alelles 
        min_diff = min(abs(len(variant[2])-len(variant[3])),abs(len(var[2])-len(var[3]))) 
        
        #set pos, start, run_end, ref, and alt_list variables
        seq = var[5]
        pos = int(var[1])
        start = int(var[4])
        run_end = start + len(seq)
        ref = var[2].upper()
        ref_end = pos + len(ref)
        alt_list = var[3].upper().split(',')
        
        
        for alt in alt_list:
            new_seq = ''
            start_pos_distance = abs(pos - start) #difference between variant and start position of variant  
              
            if pos - start >= 0: # if the starting position is within the microsatellite
                if len(ref) > len(alt): # if the reference sequence is longer than the alternate sequence
                    new_seq = seq[:start_pos_distance] + alt[:run_end - pos] + seq[start_pos_distance + len(ref):]
                else: # if the alternate sequence is longer than or equal to the reference sequence
                    new_seq = seq[:start_pos_distance] + alt[:run_end - pos] + seq[start_pos_distance + len(alt):]
         
            else: # if the starting position is before the microsatellite
                if len(alt)-len(ref) < 0:  #if the reference seq is longer than the alt
                    
                    alt_microsat_start_index = ref_end-start-min_diff # index of the start of the microsatellite in the alternate sequence
                    
                    #If the reference sequence ends before the microsatellite the microsatellite sequence is not affected by the 
                    #variant so make the new microsatellite sequence the same as the original microsatellite sequence
                    if ref_end <= start: 
                        new_seq = seq
                    else:
                        if  alt_microsat_start_index<= 0: # if the alternate allele starts before the beginning of the microsatellite
                            new_seq = seq[(ref_end - start):(ref_end - start)+len(seq)]
                        else: 
                            if ref_end>start+len(seq): # if the reference sequence ends before the end of the microsatellite
                                new_seq = alt[-alt_microsat_start_index:-(-alt_microsat_start_index+len(seq))] + seq[(ref_end - start):]
                            else: # if the reference sequence extends beyond the end of the mcirosatellite
                                new_seq = alt[-alt_microsat_start_index:] + seq[(ref_end - start):]

                else: # if the alternate sequence is longer than the refrence sequence
                    if ref_end <= start: # if the reference sequence ends before the microsatellite starts
                        new_seq = alt[len(ref):] + seq[(ref_end-start):]
                    else: # if the reference sequence does not end before the mimcrosatellite starts
                        if ref_end>start+len(seq): # if the reference sequence extends beyond the end of the microsatellite
                            new_seq = alt[-(ref_end-start+min_diff):-(ref_end-start+min_diff)+len(seq)] + seq[(ref_end - start):]
                        else: # if the reference sequence does not extend beyond the end of the microsatellite
                            new_seq = alt[-(ref_end-start+min_diff):] + seq[(ref_end - start):]                    
            # add the reconstructed sequence to the new_seq_list
            new_seq_list += [new_seq]
            
            #If there are multiple alleles on the allele list, set the sequence to the reconstructed allele and repeat the loop to
            #create a new sequence combining both of the alleles
            if c >0: 
                seq = new_seq
            c += 1

    #If there is a duplicate sequence in new_seq_list a sequence on the parent list matches the given allele, so return True. If there
    #is not duplication in new_seq_list the given allele is unique
    return len(new_seq_list) != len(set(new_seq_list)) 
    
def find_prefix_suffix(chrom,start,stop,ref_genome):
#    Find the sequences before (prefix) and after (suffix) each microsatellite using the conensus sequence and microsatellite table.
#    Input: 
#        chrom (str) - chromosome number
#        start (int) - microsatellite start position
#        stop (int) - microsatellite stop position 
#        ref_genome (Fasta object) - reference genome sequence
#    Output:
#        Prefix and suffix sequences as a pandas Series

    try:
        chrom = 'chr' + chrom[3:]
        prefix = str(ref_genome[chrom][(start-1-n*2):start+n])
        suffix = str(ref_genome[chrom][stop -n:(stop + n*2)])
        prefix = prefix.upper()
        suffix = suffix.upper()
        return pd.Series({'prefix': prefix,'suffix': suffix})
    except (ValueError, KeyError):  
        return pd.Series({'prefix': "",'suffix': ""})
    
def match_microsats_add_rows(microsatellite_table, seq_name, rf_start, rf_stop, rf_motif, rf_motif_stand,rf_prefix, rf_suffix, rf_len, rf_seq,run_name, genome_name):
#    For each microsatellite if the rf (run file) table, return the microsatellite row if it is present in the microsatellite_table.
#    
#    Input: 
#        microsatellite_table (pandas DataFrame) - master table of microsatellites
#        seq_name (str) - chromosome number
#        rf_start (int) - start index of microsatellite in run file
#        rf_stop (int) - end index of microsatellite in run file
#        rf_motif (str) - microsatellite motif in run file
#        rf_motif_stand (str) - standardised microsatellite motif in run file
#        prefix (str) - run file microsatellite prefix 
#        suffix (str) - run file microsatellite suffix
#        rf_len (int) - length of microsatellite in run file
#        rf_seq (str) - sequence of microsatellite in run file
#        run_name (str) - name of fly line being processed
#    Output:
#        seq_name (str) - chromosome number
#        rf_motif (str) - microsatellite motif in run file
#        rf_start (int) - start index of microsatellite in run file
#        rf_stop (int) - end index of microsatellite in run file
#        rf_len (int) - length of microsatellite in run file
#        prefix (str) - run file microsatellite prefix 
#        suffix (str) - run file microsatellite suffix
#        rf_seq (str) - sequence of microsatellite in run file
#        str - either 'match' or 'no_match' depending on whether a matching microsatellite was found in the master microsatellite table
#        match_ind (int) - index of matched row 

    microsatellite_table_rows = microsatellite_table.loc[abs(rf_start - microsatellite_table['SSR_Start' ]) <= sd]  #find all rows where rf_start is within sd (ex: 1500) of SSR_Start
    microsatellite_table_rows = microsatellite_table_rows.loc[microsatellite_table['Seq_Name'].str.upper() == seq_name.upper()]  #find all rows where Seq_Name matches.
    match_ind = None
    match_found = False
    
    #Go through each mcirosatellite that was pulled from the microsatellite_table and check if the rf microsatellite matches it
    for index,row in microsatellite_table_rows.iterrows():
        genomic_prefix = row[genome_name + '_prefix']
        genomic_suffix = row[genome_name + '_suffix']
        if (rf_motif_stand == row.Motif_Standardised) and (pd.isnull(genomic_prefix) != True) and (pd.isnull(genomic_suffix) != True) and \
            ((genomic_prefix[n:n*m] in rf_prefix or  genomic_prefix[:n] in rf_prefix) and (genomic_suffix[n*m:n*(m+1)] in rf_suffix or  genomic_suffix[n:n*m] in rf_suffix)):
                
            #The microsatellite is a match if Motif_Standardised matches, the prefix and suffix are not null, the segment of the microsatellite_table prefix 
            #from n to n*m (or the first n bases of the microsatellite_table prefix) is found within the rf prefix and the segment of the microsatellite_table suffix from n*m to n*(m+1) 
            #(or the segment from n to n*m) is found within the rf suffix. If the row is found in the microsatellite_table set match_found to True. If not, continue.
            match_found = True
            match_ind = row.ind

        else:
            continue
        
    #If a match was found 
    if match_found == False:

        new_row = [seq_name, rf_motif,  rf_motif_stand,rf_start, rf_stop, rf_len, rf_prefix, rf_suffix, \
            rf_seq, 'no_match', match_ind]
        return new_row

    new_row = [seq_name, rf_motif,  rf_motif_stand, rf_start, rf_stop, rf_len, rf_prefix, rf_suffix, \
        rf_seq, 'match', match_ind]
    return new_row

def find_avg(filt_cv_view, Seq_Name, SSR_Start, SSR_End):    
#    Finds the average coverage over a microsatellite region. 
#    
#    Input:
#        file_cv_view (pandas DataFrame) - portion of the coverage table 
#        Seq_Name (str) - chromosome number
#        SSR_Start (int) - start index of microsatellite
#        SSR_End (int) - end index of microsatellite
#    Output:
#        avg (float) - average coverage over the microsatellite region

    if SSR_Start == float('NaN') or SSR_End == float('NaN'):
        return 0
    if len(filt_cv_view.loc[(filt_cv_view.loc[:,('Stop')] >= int(SSR_Start)) & (filt_cv_view.loc[:,('Start')] <= int(SSR_End)),('Coverage')]) > 0:
        
        #Averages any coverage intervals that overlap the region of interest
        avg = filt_cv_view.loc[(int(SSR_Start) <= filt_cv_view.loc[:,('Stop')]) & (int(SSR_End) >= filt_cv_view.loc[:,('Start')]),('Coverage')].mean() 
        return avg
    else:
        return 0
    

def add_offspring_variants_and_match(variant_table, coverage_table, chrom,seq, par_seq_list, start, stop, motif,par_var_list, MQM):
#    For each microsatellite in the microsatellite_table uses variants within the microsatellite region to reconstruct alleles present
#    in the offspring microsatellite. Returns information about the microsatellite to add to the master microsatellite
#    table. Includes a boolean indicating whether there are alleles in the offspring microsatellite that do not match any alleles in the
#    parent microsatellite sequences.
#    
#    Input: 
#        variant_table (pandas DataFrame) - table of variants 
#        coverage_table (pandas DataFrame) - table of coverage values 
#        chrom (str) - chromosome number
#        seq (str) - microsatellite sequence
#        par_seq_list (list of strings) - list of reconstructed parent sequences for microsatellite
#        start (int) - start index of microsatellite
#        stop (int) - end index of microsatellite
#        motif (str) - microsatellite motif
#        par_var_list (list of variants) - list of variants in the parent microsatellite region
#        MQM (int) - mapping quality
#    Output:
#        all_par_match (boolean) - True if the reconstructed offspring microsatellite sequence matches any of the reconstructed parent microsatellite sequences
#        variants (list of variants) - all variants within the microsatellite region 
#        cov_avg (float) - average coverage over microsatellite region
#        bp_change_list (list of integers) - number of basepairs each variant changes the microsatellite length 
#        TYPE_list (list of strings) - type of each of the variants (snp, complex, indel, etc) 
#        genotype_list (list of strings) - genotype of each of the variants ('0/1','1/1','1/2', etc) 
#        list(set(microsat_seq_list)) (list of strings) - set of unique reconstructed sequences for the microsatellite
#        par_match_list (list of booleans) - boolean values indicating whether each reconstructed allele matches any of the reconstructed parent alleles
#        bp_change_mult_motif_list (list of integers) - list of bp changes for each allele divided by motif length, if 0 the change is a multiple of the motif 
#        AD_list (list of strings) - list of allele depth values ('1,3' indicates there is one copy of the reference allele and 2 copies of the variant allele)
#        MQM_float (list of floats) - list of MQM values for each allele
    
    #Initialize variables
    genotype_list = [] # list of the variants' genotypes    
    bp_change_list = [] # list of the number of bp and insertion or deletion changes the microsatellite length by
    par_match_list = []  # list of boolean variables indicating whether the mircrosatellite sequence matches one of the parents
    variants = ()
    TYPE_list = [] # list of boolean variables indicating whether a variant is an insertion/deletion or not.
    AD_list = []
    MQM_list = []
    cov_avg = 0.0
    seq = seq.upper()
    bp_change_mult_motif_list = []  # list indicating whether each bp change is a multiple of the motif length
    par_match = 'N/A' # indicates wether a microsatellite sequences matches one of the parent sequences
    microsat_seq_list = [seq,seq]  # use the variant to create the alternate microsatellite sequence

    #Calculate the average coverage over the mcirosatellite region
    if cov_type == 'bg':   
        cov_list = coverage_table[(start <= coverage_table.Stop) & (stop >= coverage_table.Start)]['Coverage'].tolist()
    elif cov_type == 'd':
        cov_list = coverage_table[(start <= coverage_table.POS) & (stop >= coverage_table.POS)]['Coverage'].tolist() # all the rows in the coverage table within the microsatellite region
    if cov_list != []:
        cov_list = [float(x) for x in cov_list]
        cov_avg = float(sum(cov_list)/len(cov_list)) 
    if cov_avg < cov_threshold:
        return 'I/C', variants, cov_avg,  bp_change_list, TYPE_list, genotype_list, ['I/C'], \
            par_match_list, bp_change_mult_motif_list, AD_list, MQM_list
    
    #Pull the variants out that are located between the start and stop positions of the microsatellite
    variant_rows = variant_table[(variant_table['POS'] <= stop) & ((variant_table['POS'] +variant_table['REF'].str.len()) >= start) &\
        (variant_table['#CHROM'].str.upper() == chrom.upper())][['#CHROM','POS','REF','ALT','DP','genotype','GL','TYPE', \
        'AD', 'MQM']]
    
    #If no matching variants were pulled from the offspring_vars table, set variant_list to 'no variants' and indicate that there was no 
    #variation found between the offspring microsatellite and the parent microsatellite
    if variant_rows.empty:
        par_match = seq in par_seq_list # par_match is True if the unaltered sequence is in the list of parent sequences
        par_match_list += [par_match]
        cov_list = [float('inf')]

        variants = ('no variants',)
        if False in par_match_list:
            all_par_match = False
        else:
            all_par_match = True
        if 'I/C' in par_seq_list:
            all_par_match = 'I/C'
        return all_par_match, variants, cov_avg,  bp_change_list, TYPE_list, genotype_list, \
            list(set(microsat_seq_list)),par_match_list, bp_change_mult_motif_list, AD_list, MQM_list

    
    else: # If variant_rows is not empty, pull the variant and coverage information from the rf table
        cov_list = []

            
        #Go through each variant in the microsatellite and reconstruct each alternate allele using the variant information
        for i in range(len(variant_rows)):
            current_var = variant_rows.iloc[i].tolist() # current variant stored as a list
            cov = current_var[4]
            AD_list += [str(current_var[8])]
            MQM = str(current_var[9])
            MQM_list += [MQM]
            AD_int = [int(x) for x in str(current_var[8]).split(',')]
            MQM_float = [float(x) for x in MQM.split(',')]
            pos = int(current_var[1]) # position of variant            
            cov_list += [int(cov)]
            ref = current_var[2].upper() # reference allele
            alt = current_var[3].split(',') # alternate allele
            alt = [a.upper() for a in alt]
            genotype = current_var[5]
            genotype_list += [genotype]
            
            #For each alternate allele in the variant reconstruct the microsatellite and add it to the microsat_seq_list. Also, add a boolean
            #variable to the par_match list indicating whether that allele matches one of the parents.'''
            
            microsat_seq_list += [seq]*len(alt)
                
            for j in range(len(alt)):
                
                ad = AD_int[j+1]
                mqm = MQM_float[j]
                if ad < ad_threshold:
                    variants = variants + ('I/AD',)  # insufficient allele depth
                    continue
                if mqm < mqt:
                    variants = variants + ('I/Q',)
                    continue
                variants = variants + (tuple(current_var[:3]) + (alt[j],start,seq,cov),)
                        
                bp_change_list += [len(alt[j]) - len(ref)]
                TYPE_list += [current_var[7]]

                start_pos_distance = abs(pos - start) #difference between variant and start position of variant                
                
                if pos - start >= 0: # if the starting position is within the microsatellite
                    microsat_seq_list[j+1] = microsat_seq_list[j+1][:start_pos_distance] + alt[j] + microsat_seq_list[j+1][start_pos_distance + len(ref):]

                else: # if the starting position is before the microsatellite
                    microsat_seq_list[j+1] = alt[j][start_pos_distance:] + microsat_seq_list[j+1][(len(ref) - start_pos_distance):]
    par_var_list_no_insufficients = tuple((x for x in par_var_list if x[0] != 'I'))
    variants_no_insufficients = tuple((x for x in variants if x[0] != 'I'))

    all_par_match = is_equivalent(par_var_list_no_insufficients,variants_no_insufficients)

    if (all_par_match == 'Empty Allele List'):
        if seq in par_seq_list:
            all_par_match = True
        else:
            all_par_match = False
        
    if ('I/C' in par_var_list):
        all_par_match = 'I/C'
    bp_change_mult_motif_list = [x%len(motif) for x in bp_change_list if x!=0] # keeps track of whether each variable is a multiple of the motif
    
    #Return the information about the non-matched variant that was found in the rec
    return all_par_match, variants, cov_avg,  bp_change_list, TYPE_list, genotype_list, \
        list(set(microsat_seq_list)), par_match_list, bp_change_mult_motif_list, AD_list, MQM_float
        
        
def add_parent_variants(variant_table, coverage_table, chrom,start,stop,seq, par_seq_list, par_var_list, MQM):
#    For each microsatellite in the microsatellite_table uses variants within the microsatellite region to reconstruct alleles present
#    in the parent microsatellite. Returns information about the microsatellite to add to the master microsatellite
#    table. Includes an updated list of parent variants and a list of reconstructed sequences. 
#    
#    Input:  
#        variant_table (pandas DataFrame) - table of variants 
#        coverage_table (pandas DataFrame) - table of coverage values 
#        chrom (str) - chromosome number
#        start (int) - start index of microsatellite
#        stop (int) - end index of microsatellite
#        seq (str) - microsatellite sequence
#        par_seq_list (list of strings) - list of reconstructed parent sequences for microsatellite
#        par_var_list (list of variants) - list of variants in the parent microsatellite region
#        MQM (int) - mapping quality
#    Output: variants, cov_avg, genotype_list, list(set(microsat_seq_list)), AD_list, tuple(par_var_list) + variants, MQM_float
#        variants (list of variants) - all variants within the microsatellite region 
#        cov_avg (float) - average coverage over microsatellite region
#        bp_change_list (list of integers) - number of basepairs each variant changes the microsatellite length 
#        genotype_list (list of strings) - genotype of each of the variants ('0/1','1/1','1/2', etc) 
#        list(set(microsat_seq_list)) (list of strings) - set of unique reconstructed sequences for the microsatellite
#        AD_list (list of strings) - list of allele depth values ('1,3' indicates there is one copy of the reference allele and 2 copies of the variant allele)
#        par_var_list (list of variants) - updated list of all parent variants for the microsatellites
#        MQM_float (list of floats) - list of MQM values for each allele

    cov_avg = 0.0
    
    #Pull out all rows from the cv table that are in between the start and stop positions of the microsatellite
    if cov_type == 'd':
        cov_list = coverage_table[(int(start) <= coverage_table.POS) & (int(stop) >= coverage_table.POS)]['Coverage'].tolist()
    if cov_type == 'bg':
        cov_list = coverage_table[(int(start) <= coverage_table.Stop) & (int(stop) >= coverage_table.Start)]['Coverage'].tolist()
        
    if cov_list != []:
        #average the coverages if the cov_list is not empty
        cov_list = [float(x) for x in cov_list]
        cov_avg = float(sum(cov_list)/len(cov_list))

    AD_list = [] 
    genotype_list = []
    MQM_list = []
    seq = seq.upper()
    variants = ()
    microsat_seq_list = [seq,seq]

    if cov_avg < cov_threshold:
        return ('I/C','I/C'), cov_avg, genotype_list, ['I/C'] + par_seq_list,AD_list, ('I/C',) + par_var_list,  MQM_list

    #Pull out all rows from the variant_table where the position of the variant is between the start and stop positions of the microsatellite
    variant_rows = variant_table[(variant_table['POS'] <= stop) & ((variant_table['POS']+variant_table['REF'].str.len()) >= start) & (variant_table['#CHROM'].str.upper() == \
        chrom.upper())][['#CHROM','POS','REF','ALT','DP','genotype','GL','AD','MQM']]

    if variant_rows.empty:
        cov_list = [float('inf')]
        variants = ('no variants',)
        return variants, cov_avg, genotype_list, list(set(microsat_seq_list)) + par_seq_list, AD_list, par_var_list, MQM_list 
    else:
        cov_list = []
            
        for i in range(len(variant_rows)):
            current_var = variant_rows.iloc[i].tolist() # current variant stored as a list
            cov = current_var[4]
            MQM = str(current_var[8])
            MQM_list += [MQM]
            MQM_float = [float(x) for x in MQM.split(',')]
            pos = int(current_var[1]) # position of variant
            cov_list += [int(cov)]
            ref = current_var[2].upper() # reference allele
            alt = current_var[3].split(',') # alternate allele
            alt = [a.upper() for a in alt]
            genotype = current_var[5]
            genotype_list += [genotype]
            AD_list += [str(current_var[7])]
            AD_int = [int(x) for x in str(current_var[7]).split(',')]
            microsat_seq_list += [seq]*len(alt)
                
            #For each alternate allele in the variant reconstruct the microsatellite and add it to the microsat_seq_list. Also, add a boolean
            #variable to the par_match list indicating whether that allele matches one of the parents.
            for j in range(len(alt)):

                ad = AD_int[j+1]
                mqm = MQM_float[j]
                if ad < ad_threshold:
                    variants = variants+ ('I/AD',)
                    continue
                
                if mqm < mqt:
                    variants = variants + ('I/Q',)
                    continue
                else:
                   variants = variants + (tuple(current_var[:3]) + (alt[j],start,seq,cov),)
                    
                start_pos_distance = abs(pos - start) #difference between variant and start position of variant                
                if pos - start >= 0: # if the starting position is within the microsatellite
                    microsat_seq_list[j+1] = microsat_seq_list[j+1][:start_pos_distance] + alt[j] + microsat_seq_list[j+1][start_pos_distance + len(ref):]
                else: # if the starting position is before the microsatellite
                    microsat_seq_list[j+1] = alt[j][start_pos_distance:] + microsat_seq_list[j+1][(len(ref) - start_pos_distance):]

    microsat_seq_list += par_seq_list
    return variants, cov_avg, genotype_list, list(set(microsat_seq_list)), AD_list, tuple(par_var_list) + variants, MQM_float

def calculate_matched_non_matched (microsatellite_table_dict):
#        Calculate number of matched/non-matched microsatellites for each type of variant
#        
#        Input: 
#            microsatellite_table_dict (dict) - total vs non-matched numbers of microsatellites
#        Output:
#            total_values (dict) - dictionary of summary values (total numbers)
#            percent_values (dict) - dictionary of summary values (percentages)
#            total variants (int) - total number of variants found
        
    variant_types = ['ins','del','snp','complex'] # different types of variants
    percent_type_dict = {'percent_total_': 'total_','percent_':'non_matched_'} # dictionary connecting percent type to the name of the associated total value used
    total_values = {} # dictionary to keep track of total number of variants for each variant type
    percent_values = {} # dictionary to keep track of variant percentage for each variant type
            
    for var_type in variant_types:
        
        #Calculate the total number of variants and total number of non-matched variants for each type of variant'''
#            variant_type_total = 0
        for match_type in microsatellite_table_dict:
            total_values[match_type + var_type] = microsatellite_table_dict[match_type][run_name + '_TYPE'].apply(lambda x: list(x).count(var_type)).sum()
#                variant_type_total += total_values[match_type + var_type]
            
    for match_type in microsatellite_table_dict:
        total_values[match_type + 'indels'] =  total_values[match_type  + 'ins'] + total_values[match_type  + 'del']

    total_variants = total_values['total_ins'] + total_values['total_del'] + total_values['total_snp'] + total_values['total_complex']
    
    for var_type in variant_types:
        #Calculate percentage of matched/non-matched microsatellites for each type of variant'''
        for percent_type in percent_type_dict:
            total_var_num = total_values[percent_type_dict[percent_type]+ var_type]
            percent_values[percent_type + var_type] = round((safe_div(total_var_num,total_variants))*100,2)    
    
    #Add insertions and deletions to get indel values          
    for percent_type in percent_type_dict:
        percent_values[percent_type + 'indels'] =  percent_values[percent_type + 'ins'] + percent_values[percent_type + 'del']

    return total_values,percent_values, total_variants

def calc_motif_mult (mt):
#        Calculate the percentage of changes in microsatellite length are multiples of the motif length.
#        Input:
#            mt (DataFrame) - master microsatellite table
#        Output:
#            percent_mult_motif (float) - percentage of bp changes that are multiples of the motif length
        total_bp_change = mt[run_name + '_bp_change_mult_motif'].apply(lambda x: len(list(x))).sum()
        total_num_mult_motif = mt[run_name + '_bp_change_mult_motif'].apply(lambda x: list(x).count(0)).sum()
        percent_mult_motif = round((safe_div(total_num_mult_motif,total_bp_change))*100,2)
        
        return percent_mult_motif
    
def calculate_summary_data (microsatellite_table, run_name, genotype_all, avg_cov):
#    Summarize data from the microsatellite table. Returns number and percentage of matched vs not-matched microsatellites broken down by each variant type,\
#    genotypes, total number of variatns, etc.
#    
#    Input:
#        microsatellite_table (DataFrame) - master table with information about the microsatellites and variants
#        run_name (str) - name of the fly line being analyzed
#        genotype_all (list) - number of variants with different genotypes ('0/1','1/1', etc)
#        avg_cov (float) - average coverage over all the microsatellites
#    Output:    
#        new_row (list) - list of summary data from the microsatellite table
        
    full_run_name = run_name
    run_name = run_name.split('_')[0]

    total_microsats = len(microsatellite_table.index) 
    matched_microsatellites  = microsatellite_table[(microsatellite_table[run_name + '_match'] == True) & (microsatellite_table[run_name + '_coverage_avg'] >= cov_threshold)]
    non_matched_microsatellites  = microsatellite_table[(microsatellite_table[run_name+ '_match'] == False) & (microsatellite_table[run_name + '_coverage_avg'] >= cov_threshold)]

    #Calculate the total number of matched microsatellites, non-matched microsatellites, TYPEs, bp_change, etc
    num_matched = len(matched_microsatellites.index)    
    num_not_matched = len(non_matched_microsatellites.index)
    total_matched_not_matched = num_matched + num_not_matched

    percent_matched = round((safe_div(num_matched,total_matched_not_matched))*100,2)
    percent_not_matched = round((safe_div(num_not_matched,total_matched_not_matched))*100,2)
    
    microsatellite_table_dict = {'non_matched_': non_matched_microsatellites, 'total_': microsatellite_table}
    total_values, percent_values, total_variants = calculate_matched_non_matched (microsatellite_table_dict)
    
    #Calculate percentage of motif multiple variants
    total_percent_mult_motif = calc_motif_mult(microsatellite_table)
    not_matched_percent_mult_motif = calc_motif_mult(non_matched_microsatellites)
        
    #Calculate the number of non-matched microsatellites variants of each genotpe and the average coverage for the run
    
    genotype_non_matched = calculate_genotypes(non_matched_microsatellites,run_name)
    
    new_row = [full_run_name,total_microsats, num_matched,num_not_matched,percent_matched, percent_not_matched, total_values['non_matched_ins'] + total_values['non_matched_del'], total_values['non_matched_snp'], \
        total_values['non_matched_complex'], percent_values['percent_indels'],percent_values['percent_snp'], percent_values['percent_complex'], total_percent_mult_motif, not_matched_percent_mult_motif, \
        total_values['total_ins'] + total_values['total_del'], total_values['total_snp'],total_values['total_complex'], total_variants, percent_values['percent_total_indels'], percent_values['percent_total_snp'], \
        percent_values['percent_total_complex'], avg_cov, genotype_non_matched, genotype_all, total_microsats]
    
    return new_row

def calculate_genotypes(mt,run_name):
#    Calculate the number of non-matched microsatellites variants of each genotpe and the average coverage for the run
#    Input:
#        mt (DataFrame) - microsatellite or variant table 
#        run_name (str) - name of the fly line being analyzed
#    Output:
#        genotypes (list) - list of the number of variants with different genotypes (homozygous reference, heterozygous, homozygous alternate, other)
    
    if run_name == 'variant_table':
        hom_ref = len(mt[mt.genotype == '0/0'])
        het= len(mt[mt.genotype == '0/1']) + len(mt[mt.genotype == '1/0'])
        hom_alt = len(mt[mt.genotype == '1/1'])
        other_genotypes = len(mt) - hom_ref - het - hom_alt
        genotypes = [hom_ref, het, hom_alt, other_genotypes]
        return genotypes
    else:
        hom_ref = mt[run_name + '_genotype'].apply(lambda x: x.count('0/0')).sum() # homozygous reference
        het = mt[run_name + '_genotype'].apply(lambda x: x.count('0/1')).sum() + mt[run_name + '_genotype']\
            .apply(lambda x: x.count('1/0')).sum() # heterozygous
        hom_alt = mt[run_name + '_genotype'].apply(lambda x: x.count('1/1')).sum() # homozygous alternate
        other_genotypes = mt[run_name + '_genotype'].apply(lambda x: len(list(x))).sum() -  hom_ref - het - hom_alt
        genotypes= [hom_ref, het, hom_alt,other_genotypes]
        return genotypes
    

def create_microsatellite_table (run_name, microsat_file,genome):
#    Create the microsatellite table using the file and reference genome.
#    Input:
#        run_name (str) - name of the fly line being analyzed
#        microsat_file (str) - name of the file containing the microsatellite data
#        genome (fasta object) - consensus sequence of reference genome
#    Output:
#        mt (DataFrame) - master table with information about the microsatellites and variants
    if shorten > 0:
        mt = pd.read_csv(microsat_file,sep='\t',index_col=None, nrows = shorten)
    else:
        mt = pd.read_csv(microsat_file,sep='\t',index_col=None)
    mt[['Sequence','Seq_Name']].apply(lambda x: x.fillna('',inplace=True))
    mt = mt[~mt.Motif.isnull()]
    
    #Eliminate any microsatellites that are greater than 100bp or less than 1bp long and those that do not have a valid Seq_Name
    mt = mt[(mt.Seq_Name != 'Seq_Name') & (mt['Motif'].map(len) >1) & (mt['Sequence'].map(len) <= 100)& (mt['Seq_Name'].map(len) <= 100)]
    
    #Convert the Start and End columns to numeric values and Seq_Name to categories
    mt[['SSR_Start','SSR_End']] = mt[['SSR_Start','SSR_End']].apply(lambda x: pd.to_numeric(x))
    mt.Seq_Name = mt.Seq_Name.astype('category')
    
    #Add prefix and suffix to microsatellite_table column
    mt[[run_name + '_prefix',run_name + '_suffix']] = mt.apply(lambda row: find_prefix_suffix(row['Seq_Name'], int(row['SSR_Start']), \
      int(row['SSR_End']),genome), axis=1).dropna()
    
    #Remove any rows where the prefix and suffix were not found from the table
    mt = mt[(mt[run_name + '_prefix'].map(len) >1) & (mt[run_name + '_suffix'].map(len) >1)]
    
    #Add empty new columns to the table
    mt[run_name + '_seq'] = mt['Sequence']
    mt['par_seq_list'] = [[] for x in mt['Sequence']] # list of the microsatellite sequences for each parent
    mt['par_var_list'] = [() for x in mt['Sequence']] # list of the microsatellite variants for each parent
    mt = mt.drop(['SearchMode','Mismatches','Score'],axis=1)
    return mt

def create_variant_file(variant_file):
#    Convert the variant file to a pandas table and set table columns. Truncate if shorten > 0.
#    Input:
#        variant_file (str) - name of the file containing the variant information
#    Output:
#        rf (DataFrame) - table of variant information pulled from variant_file
#        genotype (list) - number of variants with each genotype ('0/1', '1/1', etc)
    
    col_names = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','MORE_INFO']
    if shorten > 0:
        rf = pd.read_csv(variant_file,sep='\t',comment = "#", nrows = shorten, names = col_names, dtype = {'POS':int})
    else:
        rf = pd.read_csv(variant_file,sep='\t',comment = "#", names = col_names, dtype = {'POS':int})

    #Split genotype and GL columns
    rf['genotype'],rf['DP'],rf['AD'],rf['RO'],rf['QR'],rf['AO'], rf['QA'],rf['GL'] = rf.iloc[:,-1].str.split(':').str
    
    #Split the INFO column up into separate colums
    rf['TYPE'] = rf['INFO'].str.split('=').str[-1]
    rf['MQM'] = rf['INFO'].str.split(';MQM=').str[1].str.split(';').str[0]

    #Calculate the total number of variants of each genotype   
    genotype = calculate_genotypes(rf,'variant_table')
    
    return rf, genotype

def create_coverage_table (coverage_file):
#    Create coverage table.
#    Input:
#        coverage_file (str) - name of the file containing the coverage information
#    Output:
#        cv (DataFrame) - table containing coverage information
    
    if cov_type == 'bg':
        col_names = ['Seq_Name','Start','Stop','Coverage']
        column_types = {'Seq_Name':'category','Start':int, 'Stop':int,'Coverage':int}

    elif cov_type == 'd':
        col_names = ['Seq_Name','POS','Coverage']
        column_types = {'Seq_Name':'category','POS':int, 'Coverage':int}
    
    if shorten > 0:
        cv = pd.read_csv(coverage_file,sep='\t',index_col=None, nrows = shorten,header=None,dtype=column_types,names=col_names)
    else:
        cv = pd.read_csv(coverage_file,sep='\t',index_col=None,header=None,dtype=column_types,names=col_names)

    return cv

def clean_up_table_columns(mt):
#    Clean-up microsatellite table by getting ride of duplicates, extra columns, and converting to integers where necessary
#    Input:
#        mt (DataFrame) - microsatellite table
#    Output:
#        mt (DataFrame) - cleaned up microsatellite table
    
    #Drop duplicates and extra columns
    mt = mt.drop_duplicates(subset=['Seq_Name','Motif','Sequence','SSR_Start','SSR_End'])
    mt = mt.drop(['ind','match'],axis = 1)
    
    #Convert columns to integers
    mt[[run_name + '_start',run_name + '_stop']] = mt[[run_name + '_start',run_name + '_stop']].apply(lambda x: pd.to_numeric(x.fillna(-1).astype(int)))

    
    mt.loc[mt.par_seq_list.isnull(),'par_seq_list'] = mt.loc[mt.par_seq_list.isnull(),'par_seq_list'].apply(lambda x: []) # fill at empty cells in the par_seq_list column with []
    mt = mt[~mt[run_name + '_seq'].isnull()]
    return mt

def add_coverage_column(cv, microsatellite_table):
#    Add a coverage column to the mcirosatellite table
#    Input:
#        cv (DataFrame) - coverage table
#        microsatellite_table (DataFrame) - master table with information about the microsatellites and variants
#    Output:
#        microsatellite_table (DataFrame) - master table with information about the microsatellites and variants with coverage column
    
    microsatellite_table_categories = microsatellite_table.Seq_Name.unique()

    #Go through the microsatellite_table one category at a time and populate the coverage colum for each microsatellite
    for category in microsatellite_table_categories[:10]: 
        if cov_type == 'd':
            filt_cv_view = cv.loc[cv.Seq_Name == category,('Coverage','Seq_Name','POS')]
        elif cov_type == 'bg':
            filt_cv_view = cv.loc[cv.Seq_Name == category,('Coverage','Seq_Name','Start','Stop')]
        microsatellite_table.loc[microsatellite_table.Seq_Name == category,(cov_col)] = microsatellite_table.loc[microsatellite_table.Seq_Name == category].apply(lambda row: find_avg(filt_cv_view, row['Seq_Name'], \
                 row[run_name + '_start'], row[run_name + '_stop']), axis=1)  #working statement
    return microsatellite_table

#Parse Arguements
parser = argparse.ArgumentParser()                                               
parser.add_argument("--cov_type", "-covtype", type=str, required=True)
parser.add_argument("--genome_file", "-gf", type=str, required=True)
parser.add_argument("--genome_microsat_file", "-gmf", type=str, required=True)
parser.add_argument("--coverage_file", "-cf", type=str, required=True, nargs = '*') # type of coverage file (bg or d)
parser.add_argument("--variant_files", "-vcf", type=str, required=True, nargs='*')
parser.add_argument("--microsat_files", "-mf", type=str, required=True, nargs='*')
parser.add_argument("--buffer", "-b", type=int, required = True)
parser.add_argument("--shorten", "-s",type=int) # number of lines to truncate the tables to. Set to 0 if not truncating file.
parser.add_argument("--buffer_multiplier", "-m", type=int, required = True) #factor to multiply the buffer by when getting the prefix and suffix
parser.add_argument("--start_dist_match", "-sd", type=int, required = True) # how many bp apart the start distances can be to match the rows
parser.add_argument("--cov_threshold", "-ct", type=int, required = True) # only use variants with coverage at or above this threshold
parser.add_argument("--fasta_file", "-ff", type=str, required = True, nargs = '*')
parser.add_argument("--allele_depth","-ad", type=int, required = True)# number of occurances of an allele needed to count as a variant
parser.add_argument("--mq_threshold","-mq", type=int, required = True)# mappint quality threshold
parser.add_argument("--save","-save", type=bool, required = True)# save tables
args = parser.parse_args()


#Set variables using inputted parameters
n = args.buffer  # length of prefixes and suffixes
sd = args.start_dist_match # maximum distance beween start values
cov_threshold = args.cov_threshold # minimum coverage threshold
cov_type = args.cov_type # coverage type (d or bg)
shorten = args.shorten # number of lines to limit each table to
m = args.buffer_multiplier # factor to multiply buffer by when obtaining prefix and suffix
ad_threshold = args.allele_depth
mqt = args.mq_threshold # mapping quality threshold 
save = args.save # save tables or not
pd.options.mode.chained_assignment = None  # remove warning for setting values on copies 

#Read the contents of the genome and genome microsatellite files
genome = Fasta(args.genome_file)

#Create the totals tables to tally matched vs non-matched, genotpe, etc.
totals = pd.DataFrame(columns = ['Line','All Microsats','Matched Microsats','Not-matched Microsats','% Matched Microsats','% Not-matched Microsats','Not-matched indels', 'Not-matched snps', \
    'Not-matched complex', '% Not-matched indels','% Not-matched snps', '% Not-matched complex', '% Total Motif-multiple','% Not-matched Motif-multiple', 'Total indels', 'Total snps', 'Total complex',\
    'Total variants', '% Total indels','% Total snps', '% Total complex','avg_cov','genotype_non_matched','genotype_all','total_microsatelltes'])

totals_par = pd.DataFrame(columns = ['line','genotype_non_matched','genotype_all', 'avg_cov']) 

#Extract run name and genome name from file and add to par_list if necessary
drive, path = os.path.splitdrive(args.genome_microsat_file)
path, filename = os.path.split(path)
run_name = str(filename.split('_')[0])
genome_name = run_name

#Create microsatellite_table by pulling data from the genome microsatellite file.
microsatellite_table = create_microsatellite_table(genome_name, args.genome_microsat_file, genome)

del genome    
par_list = []
run_names = []

#Go through each run file and extract the variant information
for i in range(len(args.variant_files)):
    
    #Extract run name from file and add to par_list if necessary
    drive, path = os.path.splitdrive(args.variant_files[i])
    path, filename = os.path.split(path)
    run_name = str(filename.split('.')[-3])
    run_names += [run_name]
    
    print ("starting run: ",run_name, '    ', "--- %.2fs ---" % (time.time() - start_time))
    
    #Convert the variant file to a pandas table and set table columns. Truncate if shorten > 0
    rf,genotype_all = create_variant_file(args.variant_files[i])
    
    #Create coverage table
    coverage_table =  str(path.split('/')[-1]) + "_coverage"
    cov_col = coverage_table + "_avg"
    cv = create_coverage_table(args.coverage_file[i])
    
    #f the file is from a parent run that has not been processsed yet
    if (run_name[:3] == 'par' and run_name not in par_list): # THIS LINE ASSUMES THE PARENT LINE NAME STARTS WITH 'par' - CHANGE IF NECESSARY
        if args.microsat_files[i] != args.genome_microsat_file:
            
            consensus = Fasta(args.fasta_file[i])
            
            #Create microsatellite table
            mf = create_microsatellite_table(run_name, args.microsat_files[i],consensus) 
            
            #Add microsatellites to microsatellite_table if they are not already there
            microsatellite_table['ind'] = microsatellite_table.index
            new_rows = mf.apply(lambda row: match_microsats_add_rows(microsatellite_table, row['Seq_Name'], row['SSR_Start'],row['SSR_End'], row['Motif'], row['Motif_Standardised'], \
                row[run_name + '_prefix'], row[run_name + '_suffix'], row['Length'], row['Sequence'], run_name, genome_name), axis=1)
            new_rows = new_rows[~new_rows.isnull()].to_frame()
            new_rows.columns = ['Info']
            new_rows['Seq_Name'], new_rows['Motif'], new_rows['Motif_Standardised'], \
                new_rows[run_name + '_start'], new_rows[run_name + '_stop'], new_rows[run_name + '_length'], \
                new_rows[genome_name + '_prefix'], new_rows[genome_name + '_suffix'], new_rows[run_name + '_seq'], new_rows['match'], new_rows['ind']  \
                = zip(*new_rows['Info'].apply(lambda x:x[:])) # Split the Info column into separate columns
            
            #Convert rows to int and category and pull out all of the matched rows (rows where a matching microsatellite was found in the microsatellite_table)
            new_rows = new_rows.drop(['Info'],axis = 1)
            new_rows.Seq_Name = new_rows.Seq_Name.astype('category')
            new_rows_match = new_rows[new_rows.match == 'match']
                        
            #Add matched rows to the microsatellite_table if they do not already exist and clean up columns
            microsatellite_table = microsatellite_table.set_index("ind").combine_first(new_rows_match.set_index("ind")).reset_index()
            microsatellite_table = clean_up_table_columns(microsatellite_table)
            
        #Add a coverage column to the microsatellite_table (for PAR runs only - REC runs have coverage column added with add_par_match_col() )
        microsatellite_table = add_coverage_column(cv,microsatellite_table)
        par_list += [run_name] #add the run to par_list
        
        #Create new columns in the microsatellite_table
        microsatellite_table[run_name + '_coverage_avg'] = 0.0
        for column_suffix in ['_variants','_genotype','_AD','_MQM']:
            microsatellite_table[run_name + column_suffix] = [[] for x in microsatellite_table['Sequence']] 

        #Populate newly made columns using the add_offspring_variants_and_match
        microsatellite_table[run_name + '_variants'], microsatellite_table[run_name + '_coverage_avg'], microsatellite_table[run_name + '_genotype'], microsatellite_table['par_seq_list'], \
            microsatellite_table[run_name + '_AD'], microsatellite_table['par_var_list'], microsatellite_table[run_name + '_MQM'] = \
            zip(*microsatellite_table.apply(lambda row: add_parent_variants(rf, cv,row['Seq_Name'], int(row['SSR_Start']),int(row['SSR_End']), \
            row['Sequence'], row['par_seq_list'], row['par_var_list'],row[run_name + '_MQM']), axis=1).dropna())
            
        #Calculate the average coverage as well as totals for each genotype for the totals table and add values to totals table
        avg_cov = round(cv.Coverage.mean(),2)
        genotype_non_matched = calculate_genotypes(microsatellite_table,run_name)
        
        new_row = [run_name,genotype_non_matched, genotype_all, avg_cov]
        new_row_df = pd.DataFrame([new_row], columns = totals_par.columns)
        if ((totals_par.line == run_name).any()== False): # if the parent line is not already enterd into the totals table
            totals_par = totals_par.append(new_row_df) # append the new totals

        continue  # end of the parent run loop, skip the code below and proceed to the next run

#    Code for offspring runs continues here
    
    #Create new columns in the microsatellite for the variant and coverage information
    microsatellite_table[run_name + '_coverage_avg'] = 0.0
    for column_suffix in ['_variants','_genotype','_AD','_MQM']:
        microsatellite_table[run_name + column_suffix] = [[] for x in microsatellite_table['Sequence']] 
    
    #Populate newly made columns using the add_offspring_variants_and_match
    microsatellite_table = microsatellite_table[~microsatellite_table.Seq_Name.isnull()]
    microsatellite_table[run_name + '_match'], microsatellite_table[run_name + '_variants'], microsatellite_table[run_name + '_coverage_avg'], \
        microsatellite_table[run_name + '_bp_change'], microsatellite_table[run_name + '_TYPE'],microsatellite_table[run_name+ '_genotype'],microsatellite_table[run_name + '_seq_list'], \
        microsatellite_table[run_name + '_seq_par_match'],microsatellite_table[run_name + '_bp_change_mult_motif'], microsatellite_table[run_name + '_AD'],microsatellite_table[run_name + '_MQM'] \
        = zip(*microsatellite_table.apply(lambda row: add_offspring_variants_and_match(rf, cv,row['Seq_Name'],row['Sequence'], row['par_seq_list'],row['SSR_Start'],\
        row['SSR_End'],row['Motif'], row['par_var_list'],row[run_name + '_MQM']), axis=1).dropna())
    
    avg_cov = round(cv.Coverage.mean(),2)
    new_row = calculate_summary_data (microsatellite_table, run_name, genotype_all, avg_cov)
    new_row_chr2 = calculate_summary_data (microsatellite_table[(microsatellite_table.Seq_Name == 'chr2L') | (microsatellite_table.Seq_Name == 'chr2R')],run_name + '_chr2', 'NC', 'NC')
    new_row_chr3 = calculate_summary_data (microsatellite_table[(microsatellite_table.Seq_Name == 'chr3L') | (microsatellite_table.Seq_Name == 'chr3R')],run_name + '_chr3', 'NC', 'NC')
    new_row_chr4 = calculate_summary_data (microsatellite_table[(microsatellite_table.Seq_Name == 'chr4L') | (microsatellite_table.Seq_Name == 'chr4R')],run_name + '_chr4','NC', 'NC')
    new_row_chrX = calculate_summary_data (microsatellite_table[microsatellite_table.Seq_Name == 'chrX'], run_name + '_chrX', 'NC', 'NC')
    
    #Create a new row with all of the summary information from the run and append it to the totals table
    new_row_df = pd.DataFrame([new_row], columns = totals.columns)
    new_row_2_df = pd.DataFrame([new_row_chr2], columns = totals.columns)
    new_row_3_df = pd.DataFrame([new_row_chr3], columns = totals.columns)
    new_row_4_df = pd.DataFrame([new_row_chr4], columns = totals.columns)
    new_row_X_df = pd.DataFrame([new_row_chrX], columns = totals.columns)
    
    totals['Line'] = totals['Line'].astype(str)
    for t in [new_row_df,new_row_2_df,new_row_3_df,new_row_4_df,new_row_X_df]:
        totals = totals.append(t)
    
    del cv           
    del rf   

#Convert columns in the totals table to integers
totals['Matched Microsats'] = totals['Matched Microsats'].astype(int)
totals['Not-matched Microsats'] = totals['Not-matched Microsats'].astype(int)


if save:
    #Save the microsatellite_table and totals table as excel and csv files
    now = datetime.datetime.now()
    date_time=str(now.month)+ str(now.day) + str(now.hour) + str(now.minute)
    writer = pd.ExcelWriter(run_name.split('_')[0] + '_microsatellite_table_' + str(shorten) + '_' + date_time +  '.xlsx', engine='xlsxwriter')
    microsatellite_table.to_excel(writer, sheet_name='Sheet1')
    totals.to_excel(writer, sheet_name='Sheet2')
    totals_par.to_excel(writer, sheet_name='Sheet3')
    writer.save()
    #main.to_csv(run_name.split('_')[0] + '_main_table_dm6_' + str(shorten) + '_' + str(cov_threshold) + '_' + cov_type + '_'+ str(now.month) + '_' +  str(now.day)  +'_' + str(now.year)  +'_' + '.csv', sep='\t')
    #totals.to_csv(run_name.split('_')[0]+ '_' + args.totals_table[0], sep='\t')
    #totals_par.to_csv(args.totals_par_table[0], sep='\t')
    print ('saved')

print ("FINISHED" '    ', "--- %.2fs ---" % (time.time() - start_time))
print ('\n',totals_par,'\n')
print ('\n', totals)

