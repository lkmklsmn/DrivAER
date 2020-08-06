from .dca_drivaer import dca_drivaer
import tensorflow
import scanpy as sc
import pandas as pd
import os
import anndata as ad
import numpy as np
from sklearn.ensemble import RandomForestRegressor as RFR
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import numbers
import scipy


def calc_relevance(count, pheno, tf_targets, min_targets,
                   ae_type="nb-conddisp", epochs=3, early_stop=3,
                   hidden_size=(8, 2, 8), verbose=False):

    sc.pp.filter_genes(count, min_counts=1)
    gene = count.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_dca(v):
        my_counter[0] += 1
        print(f'{my_counter[0]} / {len(targets)}')

        tmp = count.copy()
        if(scipy.sparse.issparse(tmp.X)):
            tmp.X = tmp.X.toarray()
        tmp = ad.AnnData(tmp.X + 1)
        sc.pp.normalize_per_cell(tmp)
        size_factors = tmp.obs.n_counts/np.median(tmp.obs.n_counts)

        sub = count[:,v].copy()
        if(scipy.sparse.issparse(sub.X)):
            sub.X = sub.X.toarray()
        sub = ad.AnnData(sub.X + 1)
        sub.obs["size_factors"]=size_factors

        dca_drivaer(sub, mode='latent',ae_type=ae_type,epochs=epochs,
        early_stop=early_stop,hidden_size=hidden_size,verbose=verbose)
        return(sub.obsm["X_dca"])

    embed = targets.map(fun_dca)

    df_list = [pd.DataFrame(v, columns=[str(k) + '-1', str(k) + '-2']) for k, v in embed.items()]
    embed_all = None
    for df in df_list:
        embed_all = df.copy() if embed_all is None else pd.concat([embed_all, df], axis=1)

    # Random forest
    def fun_rfr(x):
        clf = RFR(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pheno)
        return rf_fit.oob_score_

    def fun_rfc(x):
        clf = RFC(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pd.factorize(pheno)[0])
        return rf_fit.oob_score_

    if isinstance(pheno[0], numbers.Number):
        rele_score = embed.map(fun_rfr)
    else:
        rele_score = embed.map(fun_rfc)

    return embed, rele_score, embed_all


def calc_relevance_pca(adata, pheno, tf_targets, min_targets):

    gene = adata.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_dca(v):
        my_counter[0] += 1
        #print(f'{my_counter[0]} / {len(targets)}')

        tmp = adata[:,v].copy()
        sc.pp.pca(tmp, n_comps= 2)

        ret = tmp.obsm['X_pca'][:,0:1]
        return(ret)

    embed = targets.map(fun_dca)

    # Random forest
    def fun_rfr(x):
        clf = RFR(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pheno)
        return rf_fit.oob_score_

    def fun_rfc(x):
        clf = RFC(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pd.factorize(pheno)[0])
        return rf_fit.oob_score_

    if isinstance(pheno[0], numbers.Number):
        rele_score = embed.map(fun_rfr)
    else:
        rele_score = embed.map(fun_rfc)

    return embed,rele_score

def calc_relevance_umap(adata, pheno, tf_targets, min_targets):

    gene = adata.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_dca(v):
        my_counter[0] += 1
        #print(f'{my_counter[0]} / {len(targets)}')

        tmp = adata[:,v].copy()
        sc.pp.pca(tmp)
        sc.pp.neighbors(tmp)
        sc.tl.umap(tmp)

        ret = tmp.obsm['X_umap']
        return(ret)

    embed = targets.map(fun_dca)

    # Random forest
    def fun_rfr(x):
        clf = RFR(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pheno)
        return rf_fit.oob_score_

    def fun_rfc(x):
        clf = RFC(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pd.factorize(pheno)[0])
        return rf_fit.oob_score_

    if isinstance(pheno[0], numbers.Number):
        rele_score = embed.map(fun_rfr)
    else:
        rele_score = embed.map(fun_rfc)

    return embed,rele_score

def calc_relevance_tsne(adata, pheno, tf_targets, min_targets):

    gene = adata.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_dca(v):
        my_counter[0] += 1
        #print(f'{my_counter[0]} / {len(targets)}')

        tmp = adata[:,v].copy()
        sc.pp.pca(tmp)
        sc.pp.neighbors(tmp)
        sc.tl.tsne(tmp)

        ret = tmp.obsm['X_tsne']
        return(ret)

    embed = targets.map(fun_dca)

    # Random forest
    def fun_rfr(x):
        clf = RFR(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pheno)
        return rf_fit.oob_score_

    def fun_rfc(x):
        clf = RFC(n_estimators=500, oob_score = True)
        rf_fit = clf.fit(X = x, y= pd.factorize(pheno)[0])
        return rf_fit.oob_score_

    if isinstance(pheno[0], numbers.Number):
        rele_score = embed.map(fun_rfr)
    else:
        rele_score = embed.map(fun_rfc)

    return embed,rele_score
