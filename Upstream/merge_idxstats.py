###A script to merge the summary output data for each sample in a project
import os
import subprocess
import sys
import pandas as pd


#make list of all files in current directory
#files = [f for f in os.listdir('./results_nt19') for if os.path.isfile(f)]
files = []
idxstats_files = []
for folder in os.listdir('./results_19nt'):
    for f in os.listdir('./results_19nt/'+folder):
        if os.path.isfile('./results_19nt/'+folder+'/'+f):
            path = './results_19nt/'+folder+'/'
            files.append(path+f)
#make list of cutadapt output files based on file naming convention
for f in files:
    if '.hsa_modbowtie_miRNA_counts.tsv' in f:
        idxstats_files.append(f)
#name output files that dataframe will be printed to
output_file = "miRbase_idxstats_summary_new.txt"

#initialize dataframe with each sample as a column and the rows for the
#counts from idxstats
df = pd.DataFrame(columns=['Sample:'])

#initialize iterating through files with idxstats summary info
for file in idxstats_files:
    with open(file, 'r') as of:
        #get the sample ID from the file name
        sample = file.split('/')[-1].split('.')[0]
        #read the idxstats into a dataframe
        of_content = pd.read_csv(of, sep='\t', header=None)
        #check if the summary dataframe is empty and if so write the genes to it
        #along with the counts for the first sample
        if df.empty:
            df['Sample:'] = of_content[0]
            df[sample] = of_content[1]
        #if genes already written to dataframe just append counts with sample ID
        else:
            df[sample] = of_content[1]
#write the summary dataframe to txt file that is tab delimited
df.to_csv(output_file,index = False, sep='\t',header=True)

