# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import subprocess


parser = argparse.ArgumentParser()                                               
parser.add_argument("--bam_file","-b", type=str, required = True,nargs=1)# bam file
parser.add_argument("--line","-l", type=str, required = True,nargs=1)# fly line
args = parser.parse_args()

pd.options.mode.chained_assignment = None  # remove warning for setting values on copies 

bam_file = args.bam_file[0]

# SET START AND STOP COORDINATE FILENAMES
microsat_start_coord_file = 'microsat_start_coord.csv'
microsat_stop_coord_file = 'microsat_stop_coord.csv'

# SET OUTPUT FILENAMES
output_file_1 = args.line[0] + 'ROI_file_1.bam'
output_file_2 = args.line[0] + 'ROI_file_2.bam'
output_file_3 = args.line[0] + 'ROI_file_3.bam'
output_file_4 = args.line[0] + '_dm6_extracted.bam'

# CREATE SAMTOOLS VIEW COMMANDS TO EXTRACT REGIONS THAT OVERLAP MICROSATELLITES
cmd1 = 'samtools view -b -h -L microsat_start_coord.csv ' + bam_file + ' > ' + output_file_1 
cmd1a = 'samtools sort ' + output_file_1 + ' > ' + output_file_2 
cmd1b = 'samtools index ' + output_file_2 
cmd2 = 'samtools view -b -h -L microsat_stop_coord.csv ' + output_file_2 + ' > ' + output_file_3
cmd2a = 'samtools sort ' + output_file_3 + ' > ' + output_file_4
cmd2b = 'samtools index ' + output_file_4

# RUN COMMANDS
subprocess.call(cmd1,shell=True)
subprocess.call(cmd1a,shell=True)
subprocess.call(cmd1b,shell=True)
subprocess.call(cmd2,shell=True)
subprocess.call(cmd2a,shell=True)
subprocess.call(cmd2b,shell=True)
