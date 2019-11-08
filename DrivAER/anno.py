import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))

def get_anno(filename,transfer):
    df = pd.read_csv(path + "/annotations/" + filename, sep='\t', header=None,
                     names=[i for i in range(2000)], low_memory=False)
    if filename == "C3.gmt":
        df = df[-df[0].str.contains("UNKNOWN")].set_index(0)
    else:
        df = df.set_index(0)
    if transfer:
        df = df.apply(lambda row: [x[0] + x[1:].lower() for x in list(row[1:].dropna())], axis=1)
    else:
        df = df.apply(lambda row: list(row[1:].dropna()), axis=1)
    return(df)



