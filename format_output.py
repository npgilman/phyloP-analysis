# mapped_data - a map containing the relationship of line in a gtf_file to all lines in a phyloP_file that map to it 
# output_path - the path the output should be written to
import pandas as pd

def format(mapped_data, phyloP_df, gtf_df, output_path):
    output_data = next(iter(output_path))
    
    header = '#chr\tstart\tend\tname\tnull_scale\talt_scale\talt_subscale\tlnlratio\tpval\n'
    df_header = header.split('\t')
    df_header[-1] = 'pval'

    new_columns = phyloP_df.columns.tolist() + df_header
    resultDf = pd.DataFrame(columns = new_columns)

    for key, arr in mapped_data.items():
        gtf_row = gtf_df.loc[key].tolist()
        for i in arr:
            phyloP_row = phyloP_df.loc[i].tolist()
            result = phyloP_row + gtf_row
            resultDf.loc[i] = result


    resultDf.to_csv(output_path, sep='\t', index=False)