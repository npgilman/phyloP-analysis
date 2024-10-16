import matplotlib.pyplot as plt
import json

speciesArr = ["Canis_lupus_familiaris", "Homo_sapiens", "Mus_musculus"]

fig, axs = plt.subplots(len(speciesArr), 22, figsize=(12*22, 5*len(speciesArr)))

speciesCnt = 0
for species in speciesArr:
    for i in range(22):
        # plt.figure(i+1)
        axs[speciesCnt, i].set_title(f"{species} chromosome {i+1}: Number of accelerated phyloP lines for each line in gtf_chr{i+1}")
        with open(f"../outputData/{species}/chr{i+1}_output.json", 'r') as file:
            data = json.load(file)
        
        vals = []
        for key, val in data.items():
            # vals.append((int(key), len(val)))
            vals.append((int(key), len(val)))
    
        x, y = zip(*vals)
        axs[speciesCnt, i].bar(x,y)
    speciesCnt = speciesCnt + 1

plt.tight_layout()
plt.savefig("accelerated_mapped_regions.png", format='png', dpi=150)
plt.show()
