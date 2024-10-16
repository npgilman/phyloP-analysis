import pandas as pd
import json

def map(chrNum, dataPath, outputPath):    
    if (chrNum == "chrX"):
        print("Aborting on chrX")
        return

    gtf_protein_data = f"./gtf_chromosmes/gtf_{chrNum}.tsv"

    proteins = pd.read_table(gtf_protein_data)
    proteinsPos = proteins[proteins.iloc[:, 6] == '+'].iloc[:, 3]
    proteinsNeg = proteins[proteins.iloc[:, 6] == '-'].iloc[:, 4]
    proteinPositions = proteinsPos.combine_first(proteinsNeg).sort_values()
    
    phyloP_data = dataPath
    
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
    
    with open(outputPath, "w") as output_file: 
        json.dump(data, output_file)
        print(f"{dataPath} mapped to {gtf_protein_data}. Written to {outputPath}")
        
