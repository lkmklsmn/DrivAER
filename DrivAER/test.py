import scanpy as sc
import pandas as pd
import os
from .rele_score import calc_relevance,rank_plot,embedding_plot
from .anno import get_anno

path = os.path.dirname(os.path.abspath(__file__))

def test_api():
    C3_mouse = get_anno(filename="C3.gmt",transfer=True)
    count1 = sc.read(path + "/data/Paul_Ery.txt", sep="\t", first_column_names=True,
                     cache=True)
    pt1 = pd.to_numeric(count1.obs_names)
    res = enrich_test(pheno = pt1, tf_targets = C3_mouse[1:5], count = count1,
                           min_targets=5, datatype = "continuous")
    rank_plot(res,save=True)
    return (res)

# enrich_plot_tf(result=res, pheno=pt1, tf_name="CMYB_01", datatype="continuous", save=False)
# enrich_plot_marker(result=res1, tf_name="YTATTTTNR_MEF2_02", gene="Zfpm1", count=count1, save=False)
