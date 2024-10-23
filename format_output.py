# mapped_data - a map containing the relationship of line in a gtf_file to all lines in a phyloP_file that map to it 
# output_path - the path the output should be written to
import pandas as pd

def format(mapped_data, phyloP_df, gtf_df, output_path):
    output_data = next(iter(output_path))

    # gtf_header = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
    new_columns = phyloP_df.columns.tolist() + gtf_df.columns.tolist()
    resultDf = pd.DataFrame(columns = new_columns)

    gtf_df.iloc[:,8] = gtf_df.iloc[:,8].apply(lambda x : x.split(';')).apply(lambda x : [y.split() for y in x])
    gtf_df.iloc[:,8] = gtf_df.iloc[:,8].apply(lambda data : [item for item in data if item and (item[0] == 'gene_id' or item[0] == 'gene_synonym')])

    for key, arr in mapped_data.items():
        gtf_row = gtf_df.loc[key].tolist()
        for i in arr:
            phyloP_row = phyloP_df.loc[i].tolist()
            result = phyloP_row + gtf_row
            resultDf.loc[i] = result

    resultDf = resultDf.drop(columns=['5', '7'])

    resultDf.to_csv(output_path, sep='\t', index=False)