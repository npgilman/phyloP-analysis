#Resources: https://bluegenes.github.io/snakemake-via-slurm/
#https://hackmd.io/@bluegenes/BJPrrj7WB

## Imports
import csv
import pandas as pd
import json

## Config file and directories
workdir: "/orange/kgraim/data/phyloP/Accelerated_regions"
log_dir = "logs"

species_list = "/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/ARs/species_list.txt"

#Read in accession
#giada runs: 5,18,21,22,X
#TREE = []
#mammals= [*range(1, 240, 1)]
# # define way to prepend the mammals, Source: https://www.geeksforgeeks.org/python-insert-the-string-at-the-beginning-of-all-items-in-a-list/
#def prepend(list, str):
#    # Using format()
#    str += '{0}'
#    list = [str.format(i) for i in list]
#    return(list)
#prefix = "fullTreeAnc"
#TREE = prepend(mammals, prefix)
#TREE.append("fullTreeAnc110point5")
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
        expand(["/orange/kgraim/data/phyloP/Accelerated_regions" + "/" + "{mammal}" + "/" + "chr" + "{chrom_num}" + "_output.tsv"], mammal=MAMMALS, chrom_num=CHROM_NUM)

## Pipeline
rule run_phyloP:
    log:
        stdout=log_dir+ "/" + "{mammal}" + "_" + "{chrom_num}.stdout",
        stderr=log_dir + "/" + "{mammal}" + "_" + "{chrom_num}.stderr"

    input:
        dataPath = expand(["/orange/kgraim/data/phyloP/Input_files" + "/" + "{mammal}" + "/" + "Enc_full.chr" + "{chrom_num}" + ".bed.acc.unfilt.phyloP"], mammal=MAMMALS, chrom_num=CHROM_NUM)
	gtf_protein_data = expand(["/blue/kgraim/ngilman/phyloP/gtf_chromosomes/gtf_chr" + "{chrom_num}" + ".tsv"], chrom_num=CHROM_NUM)

    params:
	dir = "/orange/kgraim/data/phyloP/Accelerated_regions" + "/" + "{mammal}"
        mammal = "{mammal}",
        chrom = "{chrom_num}",
        #in_fna = "/blue/james.cahill/share/241-mammalian-2020v2.chr21.maf.gz.all.fa",
        #in_mod = "/blue/james.cahill/gi.padovani/zoonomia_ARs" +"/"+  "200m_model_242.mod"

    output:
        outputPath = "/orange/kgraim/data/phyloP/Accelerated_regions" + "/" + "{mammal}" + "/" + "chr" + "{chrom_num}" + "_output.tsv"

    run:
        shell("mkdir -p {params.dir}")
	if (chrNum == "chrX"):
	        print("Aborting on chrX")
        return

        proteins = pd.read_table({input.gtf_protein_data})
        proteinsPos = proteins[proteins.iloc[:, 6] == '+'].iloc[:, 3]
        proteinsNeg = proteins[proteins.iloc[:, 6] == '-'].iloc[:, 4]
        proteinPositions = proteinsPos.combine_first(proteinsNeg).sort_values()
    
        phyloP_data = {input.dataPath}
    
        phyloP = pd.read_table(phyloP_data)
        phyloP["name"] = phyloP["name"].apply(lambda x: x.split(":")[1].split("-"))
        phyloP.rename(columns={"name" : "range"}, inplace=True)
    
        currentPosition = 0
        nextPosition = 0
    
        data = {}
        phyloPLineCounter = 0
        for positionIndex in range(proteinPositions.count()):
            currentPosition = proteinPositions[proteinPositions.index[positionIndex]]
            data[int(proteinPositions.index[positionIndex])] = []
    
            if (positionIndex + 1 < proteinPositions.count()):
                nextPosition = proteinPositions[proteinPositions.index[positionIndex + 1]]
    
            if (phyloPLineCounter >= phyloP["range"].count()):
                    print(f"{proteinPositions.index[positionIndex]} after all valid phyloP maps")
                    break
            rangeStart = int(phyloP.loc[phyloPLineCounter, :]["range"][0])
            rangeEnd = int(phyloP.loc[phyloPLineCounter, :]["range"][1])
    
            while (min(abs(currentPosition - rangeStart), abs(currentPosition - rangeEnd)) < min(abs(nextPosition - rangeStart), abs(nextPosition - rangeEnd))):
                # print(f"appending {phyloPLineCounter} to {proteinPositions.index[positionIndex]}")
                if (phyloP.loc[phyloPLineCounter, :]["lnlratio"] > 0):
                    data[proteinPositions.index[positionIndex]].append(phyloPLineCounter)
                phyloPLineCounter = phyloPLineCounter + 1
                if (phyloPLineCounter >= phyloP["range"].count()):
                    break
                rangeStart = int(phyloP.loc[phyloPLineCounter, :]["range"][0])
                rangeEnd = int(phyloP.loc[phyloPLineCounter, :]["range"][1])
    
        with open({output.outputPath}, "w") as output_file:
            json.dump(data, output_file)
        print(f"{input.dataPath} mapped to {input.gtf_protein_data}. Written to {output.outputPath}")

