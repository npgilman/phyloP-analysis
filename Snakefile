#Resources: https://bluegenes.github.io/snakemake-via-slurm/
#https://hackmd.io/@bluegenes/BJPrrj7WB

## Imports
# python packages
import csv
import pandas as pd
import json
import os
from pathlib import Path

# Files
import map_regions
import format_output

## Config file and directories
workdir: "/orange/kgraim/data/phyloP/Accelerated_regions"
log_dir = "logs"

species_list = "/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/ARs/species_list.txt"

CHROM_NUM = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
MAMMALS = []
with open(species_list, "r") as acc_file:
    csv_reader = csv.reader(acc_file)
    #next(csv_reader, None)
    for row in csv_reader:
        if (len(row) > 0):
            MAMMALS.append(row[0])


## All rule

rule all:
    input:
        expand(["/orange/kgraim/data/phyloP/Accelerated_regions/{mammal}/chr{chrom_num}_output.tsv"], mammal=["Mus_musculus"], chrom_num=CHROM_NUM)
        # dataPath = expand(["/orange/kgraim/data/phyloP/Input_files" + "/" + "{mammal}" + "/" + "Enc_full.chr" + "{chrom_num}" + ".bed.acc.unfilt.phyloP"], mammal=["Mus_musculus"], chrom_num=CHROM_NUM)

rule run_hypothesisCorrection:
    log:
        stdout=log_dir + "/correction_{mammal}_{chrom_num}.stdout",
        stderr=log_dir + "/correction_{mammal}_{chrom_num}.stderr"

    input:
        dataPath = "/orange/kgraim/data/phyloP/Input_files/cat_runs/{mammal}.phyloP",
    
    output: 
        outputPath = "/orange/kgraim/data/phyloP/Input_files/corrected_runs/{mammal}_corrected.tsv"

    params:
        mammal = "{mammal}"

    run:
        RScript <filename> {params.mammal}


## Pipeline
rule run_phyloP:
    log:
        stdout=log_dir + "/mapping_{mammal}_{chrom_num}.stdout",
        stderr=log_dir + "/mapping_{mammal}_{chrom_num}.stderr"


    # input needs to accept a single file
    input:
        dataPath = "/orange/kgraim/data/phyloP/Input_files/{mammal}/Enc_full.chr{chrom_num}.bed.acc.unfilt.phyloP",
        gtf_protein_data = "/blue/kgraim/ngilman/phyloP/gtf_chromosomes/gtf_chr{chrom_num}.tsv"

    params:
        dir = "/orange/kgraim/data/phyloP/Accelerated_regions/{mammal}",
        mammal = "{mammal}",
        chrom = "{chrom_num}"

    output:
        outputPath = "/orange/kgraim/data/phyloP/Accelerated_regions/{mammal}/chr{chrom_num}_output.tsv"

    run:
        shell("mkdir -p {params.dir}")
        data = map_regions.map({input.dataPath}, {input.gtf_protein_data}, {output.outputPath})
        # print(data)

        # current the separate format function doesn't work because it 
        # needs access to the phyloP dataframe in the map function.
        #
        # current workaround is putting the below function into the map function at the end
        # format_output.format(data, {output.outputPath})
        
	    

