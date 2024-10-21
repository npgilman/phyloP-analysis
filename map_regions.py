import pandas as pd
import json

# phyloP_path - path to a phyloP file that is read in and mapped to given gtf file
# gtf_path - path to a gtf_file that is used to map the phyloP file
def map(phyloP_path, gtf_path, output_path):    

    phyloP_path = (next(iter(phyloP_path)))
    gtf_path = (next(iter(gtf_path)))
    output_path = next(iter(output_path))


    proteins = pd.read_table(gtf_path)
    proteinsPos = proteins[proteins.iloc[:, 6] == '+'].iloc[:, 3]
    proteinsNeg = proteins[proteins.iloc[:, 6] == '-'].iloc[:, 4]
    proteinPositions = proteinsPos.combine_first(proteinsNeg).sort_values()
    
    phyloP = pd.read_table(phyloP_path)
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
    
    # with open("dump.txt", "w") as output_file: 
    #     json.dump(data, output_file)
    #     print(f"Map dumped.")


    new_columns = phyloP.columns.tolist() + ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    resultDf = pd.DataFrame(columns = new_columns)

    for key, arr in data.items():
        gtf_row = proteins.loc[key].tolist()
        for i in arr:
            phyloP_row = phyloP.loc[i].tolist()
            result = phyloP_row + gtf_row
            resultDf.loc[i] = result


    resultDf.to_csv(output_path, sep='\t', index=False)

    return data
        
