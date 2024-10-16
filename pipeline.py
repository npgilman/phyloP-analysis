import os
import map_regions

data_dir = f"../phyloPData"
species = next(os.walk(data_dir))[1]

output_dir = f"../outputData"

species_map = {}
for s in species:
    for path in next(os.walk(data_dir + "/" + s))[2]:
        chrNum = (path.split('.')[1])
        
        outputDir = (output_dir + "/" + s)
        os.makedirs(outputDir, exist_ok=True)
        
        inputData = (f"../phyloPData/{s}/Enc_full.{chrNum}.bed.acc.unfilt.phyloP")
        outputData = (f"../outputData/{s}/{chrNum}_output.json")

        map_regions.map(chrNum, inputData, outputData)
        
