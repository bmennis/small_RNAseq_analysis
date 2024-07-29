###A script to merge the summary output data for each sample in a project
import os
import subprocess
import sys
import pandas as pd


#make list of all files in current directory
files = [f for f in os.listdir('.') if os.path.isfile(f)]
cutadapt_files = []
#make list of cutadapt output files based on file naming convention
for f in files:
    if 'output-1_smallRNAseq_pipeline_' in f and '.txt' in f:
        cutadapt_files.append(f)

#name output files that dataframe will be printed to
output_file = "cutadapt_summary_new.txt"
#establish the rows containing the specific reads per sample
df_rows = [
"Total reads processed:",
"Reads with adapters:",
"Reads that were too short:",
"Reads that were too long:",
"Reads written (passing filters):"
]
#initialize dataframe with each sample as a column and the rows for the
#summarized reads for those respective samples
df = pd.DataFrame(df_rows, columns=['Sample:'])

#initialize iterating through files with cutadapt summary info
for file in cutadapt_files:
    with open(file, 'r') as of:
        of_content = of.readlines()
        #get the sample ID from the command line, specifically the output folder for that sample
        sample = of_content[4].split("results_19nt/")[-1].split('/')[0]
        #get the summary cutadapt information for the single sample
        summary = of_content[10:15]
        summary_data = []
        #iterate through each line to get just the reads numbers and not the row headers
        for line in summary:
            line = line.strip().split(' ')[-2:]
            new_line = ' '.join(line).strip()
            summary_data.append(new_line)
        #add the new column to the dataframe with the sample id as the column header
        df[sample] = summary_data
df.to_csv(output_file,index = False, sep='\t',header=True)
