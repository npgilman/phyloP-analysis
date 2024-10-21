# mapped_data - a map containing the relationship of line in a gtf_file to all lines in a phyloP_file that map to it 
# output_path - the path the output should be written to
def format(mapped_data, output_path):
    output_data = next(iter(output_path))

    new_columns = phyloP.columns.tolist() + ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    resultDf = pd.DataFrame(columns = new_columns)

    for key, arr in data.items():
        gtf_row = proteins.loc[key].tolist()
        for i in arr:
            phyloP_row = phyloP.loc[i].tolist()
            result = phyloP_row + gtf_row
            resultDf.loc[i] = result


    resultDf.to_csv(output_path, sep='\t', index=False)