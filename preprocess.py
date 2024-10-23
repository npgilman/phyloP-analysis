import pandas as pd
from io import StringIO

def preprocess(input_phyloP):
    gtf_folder = './gtf_chromosomes'
    with open(input_phyloP, 'r') as file:
        content = file.read()
    
    header = '#chr\tstart\tend\tname\tnull_scale\talt_scale\talt_subscale\tlnlratio\tpval\n'
    df_header = header.split('\t')
    df_header[-1] = 'pval'
    sections = content.strip().split(header)
    sections = sections[1:]

    dataframes = []
    for section in sections:
        df = pd.read_csv(StringIO(section), sep='\t', header=None)
        df.columns = df_header
        dataframes.append(df)

    chr_dfs = {}
    for df in dataframes:
        if (df.iloc[0,0] == 'chrX'):
            continue
        proteins = pd.read_csv(f"{gtf_folder}/gtf_{df.iloc[0,0]}.tsv", sep='\t')
        # print(proteins.head())
        chr_dfs[df.iloc[0,0]] = [df, proteins]

    return chr_dfs