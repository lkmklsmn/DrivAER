import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))

def get_anno(filename,filetype,conv_mouse):
    if filetype=="gmt":
        df = pd.read_csv(path + "/annotations/" + filename, sep='\t', header=None,
                         names=[i for i in range(2000)], low_memory=False)
        if filename == "C3.gmt":
            df = df[-df[0].str.contains("UNKNOWN")].set_index(0)
        else:
            df = df.set_index(0)
        if conv_mouse:
            df = df.apply(lambda row: [x[0] + x[1:].lower() for x in list(row[1:].dropna())], axis=1)
        else:
            df = df.apply(lambda row: list(row[1:].dropna()), axis=1)
        return(df)
    if filetype=="tsv":
        df = pd.read_csv(path + "/annotations/" + filename, sep="\t",
                         names=["TF", "Target", "Type", "Source"])
        if conv_mouse:
            df['Target'] = df['Target'].map(lambda x:x[0] + x[1:].lower())
        df = df.groupby("TF")['Target'].unique()
        return(df)


