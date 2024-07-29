### load modules needed
import os
import subprocess
import sys

sample = []
sample_full = []
batchfile = "miRNA_batch_run.sh"
temp1= """#!/bin/bash -l
#SBATCH --job-name=star
"""

temp2= """

#SBATCH --ntasks=18
#SBATCH --time=99:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH -p compute
#SBATCH --mail-type=ALL
####SBATCH mail-user=bennis@pennstatehealth.psu.edu
#SBATCH --output=output-1_smallRNAseq_pipeline_%J.txt       # Name of output file
#SBATCH --error=error-1_smallRNAseq_pipeline_err%J.txt      # Name of error msg file


echo Running on 'hostname'
echo Clean and Align is started 'date'
echo "starting"

module load umi_tools
module load cutadapt
module load bowtie
module load samtools
module load bedtools

miRNA_REF="/gpfs/Cores/GSC/share/refGenomes/miRNA/hsa_mature.bowtie"
genome_REF="/gpfs/Cores/GSC/share/refGenomes/human/ENCODE_miRNA/GRCh38_BOWTIE_Index/GRCh38_bowtie"
tag_bed="/gpfs/Cores/GSC/share/refGenomes/miRNA/hsa.mature.miRNA.bed"
refpre="hsa"
"""
### TODO:CHANGE THIS LINE BASED ON PROJECT --- EXAMPLE
project=sys.argv[1]

fastqDir="/gpfs/Cores/GSC/share/03_Analysis/" + project + "/raw_data"

resultDir="/gpfs/Cores/GSC/share/03_Analysis/" + project + "/results_19nt"

subprocess.run("mkdir -p "+resultDir, shell=True)

## Read all FASTQ files from folder
for file in os.listdir(fastqDir):
    ### different file name format, different syntax below MAY GET STUCK HERE
    file1 = file[:-16]
    fileFull = file

    if file1 in sample:
        continue
    else:
        sample_full.append(fileFull)
        sample.append(file1)
        subprocess.run("mkdir -p "+resultDir+"/"+file1+"/", shell=True)

print(sample)


####### Prepare shell scripts for each sample which will be used to seq cleaning and seq alignment
with open(batchfile,"w") as outfile:
    outfile.write("#!/usr/bin/env bash\n")
    i = 0
    for file in sample:
        file1 = file+"_R1_001.fastq.gz"
        #file2 = fastqDir+file+"_R2_001.fastq.gz"
        resultDirFileExt = resultDir + "/" + file + "/" + file
        resultDirFileExtPre = resultDirFileExt + "/" + file1
        pfile = "p1_" + file +"_miRNA_full.sh"
        ofile = "p1_" + file +"_miRNA_full.out"
        outfile.write("sbatch "+pfile+ "\n")
        with open(pfile,"w") as poutfile:
            print("###### Processing sample" + file + "#######")
            poutfile.write(temp1)
            poutfile.write("\n")
            #SBATCH -o p2.1_star_hum.out
            poutfile.write("#SBATCH -o "+ ofile + "\n")
            poutfile.write(temp2 + "\n")
            ####### Small RNAseq steps are as follows

            #umi extraction - discard AACTGTAGGCACCATCAAT on the read with a maximum of 2 mismatches; capture and keep the next 12 characters as the UMI, discard anything after the UMI tag
            #cmd1 = "umi_tools extract --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)' -L " + resultDirFileExt + "1.umi_extract.12.log.file -I " + fastqDir + "/" + sample_full[i]
+ " -S "+resultDirFileExtPre+"_umi_R1_001.fastq.gz"

#Did not run umi_tools for trimming

            #retain reads that are between 19 and 55 bp long
            cmd2 = "cutadapt -m 19 -M 55 -a AACTGTAGGCACCATCAAT -o " + resultDirFileExt+"_cutadapt_R1_001.fastq.gz "+fastqDir+"/"+file1+" 2>"+resultDirFileExt+"_p1."+file+".${refpre}.cutadapt.log.file"

            #Align reads to the known mature miRNA sequences from miRBase using Bowtie. Using Dr. Hick's relaxed parameters - allow 2 mismatches, report up to 1 alignment per read, keep unaligned reads
            cmd3 = "bowtie -p12 -v2 --best -n 3 --strata -a -k1 --un "+resultDirFileExt+"_modbowtie_R1_001.unaligned.fastq ${miRNA_REF} "+resultDirFileExt+"_cutadapt_R1_001.fastq.gz -S "+resultDirFileExt+".${refpre}.modbowtie_mature.miRNA.sam 2>
  "+resultDirFileExt+"_p2."+file+".${refpre}.miRBase_miRNA.bowtie.log"

            #Sort the SAM file and convert to BAM
            cmd4 = "samtools sort -m 4G -o "+resultDirFileExt+".${refpre}.modbowtie_mature.miRNA.sorted.bam -O bam --threads 12 "+resultDirFileExt+".${refpre}.modbowtie_mature.miRNA.sam"

            #Index the BAM file
            cmd5 = "samtools index "+resultDirFileExt+".${refpre}.modbowtie_mature.miRNA.sorted.bam"

            #Dedeuplicate the BAM file where the reads are grouped according to the exact UMI
#            cmd6 = "umi_tools dedup –method=unique -L "+resultDirFileExtPre+".${refpre}_modbowtie_miRNA_dedup.log.txt --timeit="+resultDirFileExtPre+".${refpre}_modbowtie_dedup.time.log --output-stats="+resultDirFileExtPre+".${refpre}_modbowtie
_miRNA_dedup_STATS.txt --stdin="+resultDirFileExtPre+".${refpre}.modbowtie_mature.miRNA.sorted.bam --stdout="+resultDirFileExtPre+".${refpre}.modbowtie_miRNA.dedup.bam"

            #Index the BAM file
#            cmd7 = "samtools index "+resultDirFileExt+".${refpre}.modbowtie_miRNA.dedup.bam"

#Did not run deduplicate with umi_tools because of no UMI used

            #Record the count of each miRNA in every sample
            cmd8 = "samtools idxstats "+resultDirFileExt+".${refpre}.modbowtie_mature.miRNA.sorted.bam | cut -f 1,3 > "+resultDirFileExt+".${refpre}_modbowtie_miRNA_counts.tsv"

            #Align the unaligned reads to the Human reference genome using the same parameters as above
            cmd9 = "bowtie -p12 -v2 --best -n 3 --strata -a -k1 --un "+resultDirFileExt+".${refpre}_modbowtie_genome_unaligned.fastq ${genome_REF} "+resultDirFileExt+"_modbowtie_R1_001.unaligned.fastq  -S "+resultDirFileExt+".${refpre}.modbowtie_genome.sam 2> "+resultDirFileExt+"_p4."+file+".${refpre}_genome_bowtie.log.txt"

            #Sort the SAM output file and write to a BAM file
            cmd10 = "samtools sort -m 4G -o "+resultDirFileExt+".${refpre}.modbowtie_genome.sorted.bam -O bam --threads 6 "+resultDirFileExt+".${refpre}.modbowtie_genome.sam"

            #Index the BAM file
            cmd11 = "samtools index "+resultDirFileExt+".${refpre}.modbowtie_genome.sorted.bam"

            #Deduplicate the BAM file where the reads are grouped according to the exact UMI
            #cmd12 = "umi_tools dedup –method=unique -L "+resultDirFileExtPre+".${refpre}.modbowtie_genome_dedup.log.txt --output-stats="+resultDirFileExtPre+".${refpre}.modbowtie_genome_dedup_STATS.txt --stdin="+resultDirFileExtPre+".${refpre}.modbowtie_genome.sorted.bam --stdout="+resultDirFileExtPre+".${refpre}.modbowtie_genome.dedup.bam"

            #Index the BAM file
            #cmd13 = "samtools index "+resultDirFileExtPre+".${refpre}.modbowtie_genome.sorted.bam"

#Did not run deduplicate with umi_tools because of no UMI used

            #Annotate only those regions in the deduplicated BAM file which fall within the coordinates provided by miRBase in the form of a BED file and generate a tagged BAM file using tagBAM from Bedtools software
            cmd14 = "tagBam -i "+resultDirFileExt+".${refpre}.modbowtie_genome.sorted.bam -files ${tag_bed} -names -tag XQ > "+resultDirFileExt+".${refpre}_modbowtie_genome_tagged.bam"
            #Count the miRNAs from this second step using the appropriate "get_counts.sh" script

            # --readFilesCommand zcat
            poutfile.write("\n")
            #os.system(cmd1)
#            poutfile.write(cmd1+"\n")
            poutfile.write(cmd2+"\n")
            poutfile.write(cmd3+"\n")
            poutfile.write(cmd4+"\n")
            poutfile.write(cmd5+"\n")
#            poutfile.write(cmd6+"\n")
#            poutfile.write(cmd7+"\n")
            poutfile.write(cmd8+"\n")
            poutfile.write(cmd9+"\n")
            poutfile.write(cmd10+"\n")
            poutfile.write(cmd11+"\n")
#            poutfile.write(cmd12+"\n")
#            poutfile.write(cmd13+"\n")
            poutfile.write(cmd14+"\n")
            i = i + 1
