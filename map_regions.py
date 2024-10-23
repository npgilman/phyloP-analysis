import pandas as pd

# phyloP_path - a dataframe that stores the gtf info for a chromosome
# gtf_path - path to a gtf_file that is used to map the phyloP file
def map(phyloP, proteins):   
    proteinsPos = proteins[proteins.iloc[:, 6] == '+'].iloc[:, 3]
    proteinsNeg = proteins[proteins.iloc[:, 6] == '-'].iloc[:, 4]
    proteinPositions = proteinsPos.combine_first(proteinsNeg).sort_values()
    
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

    return data
        
