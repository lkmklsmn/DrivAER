from .dca_drivaer import dca_drivaer
import tensorflow as tf
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
import matplotlib.pyplot as plt
import seaborn as sns


def calc_relevance(count, pheno, tf_targets, min_targets,
                   ae_type="nb-conddisp", epochs=50, early_stop=3,
                   hidden_size=(8, 2, 8), verbose=False):

    sc.pp.filter_genes(count, min_counts=1)
    gene = count.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]
    # "targets" is a pandas series, with index representing TF names. Each row is a list of target genes.

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

        tf.keras.backend.clear_session()

        return(sub.obsm["X_dca"])

    embed = [fun_dca(x) for x in targets]

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
        rele_score = [fun_rfr(x) for x in embed]
    else:
        rele_score = [fun_rfc(x) for x in embed]
        
    names = targets.index.tolist()

    return embed, rele_score, names
    # return embed, rele_score


def calc_relevance_pca(adata, pheno, tf_targets, min_targets):

    gene = adata.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_pca(v):
        my_counter[0] += 1
        #print(f'{my_counter[0]} / {len(targets)}')

        tmp = adata[:,v].copy()
        sc.pp.pca(tmp, n_comps= 2)

        ret = tmp.obsm['X_pca'][:,0:1]
        return(ret)

    embed = [fun_pca(x) for x in targets]

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
        rele_score = [fun_rfr(x) for x in embed]
    else:
        rele_score = [fun_rfc(x) for x in embed] 
    
    names = targets.index.tolist()
    return embed, rele_score, names
    # return embed, rele_score


def calc_relevance_umap(adata, pheno, tf_targets, min_targets):

    gene = adata.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_umap(v):
        my_counter[0] += 1
        #print(f'{my_counter[0]} / {len(targets)}')

        tmp = adata[:,v].copy()
        sc.pp.pca(tmp, n_comps= 10)
        sc.pp.neighbors(tmp)
        sc.tl.umap(tmp)

        ret = tmp.obsm['X_umap']
        return(ret)

    embed = [fun_umap(x) for x in targets]

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
        rele_score = [fun_rfr(x) for x in embed]
    else:
        rele_score = [fun_rfc(x) for x in embed]

    names = targets.index.tolist()
    return embed, rele_score, names
    # return embed, rele_score


def calc_relevance_tsne(adata, pheno, tf_targets, min_targets):

    gene = adata.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]

    my_counter = [0]

    def fun_tsne(v):
        my_counter[0] += 1
        #print(f'{my_counter[0]} / {len(targets)}')

        tmp = adata[:,v].copy()
        sc.pp.pca(tmp, n_comps= 10)
        sc.pp.neighbors(tmp)
        sc.tl.tsne(tmp)

        ret = tmp.obsm['X_tsne']
        return(ret)

    embed = [fun_tsne(x) for x in targets]

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
        rele_score = [fun_rfr(x) for x in embed]
    else:
        rele_score = [fun_rfc(x) for x in embed]

    names = targets.index.tolist()
    return embed, rele_score, names
    # return embed, rele_score


def plot_random(original_score, random_scores):

  sns.distplot(random_scores, hist=False, rug=True)
  plt.axvline(original_score, 0, 2,  color = 'red')
  plt.title('Distribution random genes')


def compare_to_random(count, pheno, geneset,
                      ae_type = 'nb-conddisp', min_targets = 10, epochs = 50, early_stop = 3,
                      num_permutations = 10, plot = True):

  res = calc_relevance(count = count,
                              pheno = pheno,
                              ae_type = ae_type,
                              tf_targets = geneset,
                              min_targets = min_targets,
                              epochs=epochs,
                              early_stop = early_stop)

  original_score = res[1][0]

  genesets = []
  for x in range(0, num_permutations):
    genesets.append(list(np.random.choice(list(count.var_names), len(geneset[0]))))
  random_genesets = pd.Series(genesets)

  random = calc_relevance(count = count,
                              pheno = pheno,
                              ae_type = ae_type,
                              tf_targets = random_genesets,
                              min_targets = min_targets,
                              epochs=epochs,
                              early_stop = early_stop)

  if plot:
    plot_random(original_score, list(random[1]))

  return original_score, list(random[1])
