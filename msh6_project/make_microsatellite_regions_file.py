# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import subprocess

parser = argparse.ArgumentParser()  
parser.add_argument("--td_file", "-td", type =str, required=True, nargs=1) # input td file                                             
args = parser.parse_args()

pd.options.mode.chained_assignment = None  # remove warning for setting values on copies 

# CLEAN UP .TD FILE TO REMOVE EXTRANEOUS MICROSATELLITES
awk_cmd = "awk \'{if ((substr($1,1,3)==\"chr\") && (int($4)>=0) && (int($5)>=0)) printf (\"%s\t%i\t%i \\n\",$1,$4,$5)}\' " + args.td_file[0] + " > microsats_clean.td" 
subprocess.call(awk_cmd,shell=True)

# IMPORT MICROSATS_CLEAN.TD FILE AS A PANDAS TABLE 
columns=['Seq_Name','SSR_Start', 'SSR_End']
main_table = pd.read_csv('microsats_clean.td',sep='\t',index_col=None,header=None,names=columns,  dtype={'SSR_End':'int'})

# CREATE TABLE OF START COORDINATES
microsat_start_coord =  main_table[['Seq_Name','SSR_Start','SSR_Start']]
microsat_start_coord.columns = ['chrom', 'start_1','start_2']
microsat_start_coord['start_1'] = microsat_start_coord['start_1']-10
microsat_start_coord['start_2'] = microsat_start_coord['start_2']-8
microsat_start_coord['start_1'][microsat_start_coord.start_1 <0] = 0 # make sure table does not contain negative coordinates
microsat_start_coord['start_2'][microsat_start_coord.start_2 <0] = 0

# CREATE TABLE OF STOP COORDINATES
microsat_stop_coord =  main_table[['Seq_Name','SSR_End','SSR_End']]
microsat_stop_coord.columns = ['chrom', 'stop_1','stop_2']

microsat_stop_coord['stop_1'] = microsat_stop_coord['stop_1']+8
microsat_stop_coord['stop_2'] = microsat_stop_coord['stop_2']+10

# SAVE START AND STOP COORDINATE FILES
microsat_start_coord.to_csv('microsat_start_coord.csv',sep="\t",index=None,header=None)
microsat_stop_coord.to_csv('microsat_stop_coord.csv',sep="\t",index=None,header=None)
