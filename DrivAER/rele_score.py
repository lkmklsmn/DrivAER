from dca.api import dca
import tensorflow
import scanpy as sc
import pandas as pd
import os
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor as RFR
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

my_counter = 0
path = os.path.dirname(os.path.abspath(__file__))

def calc_relevance(count, pheno, datatype, tf_targets, min_targets,
                   ae_type="nb-conddisp", epochs=3, early_stop=3,
                   hidden_size=(8, 2, 8), verbose=False):
    gene = count.var_names.tolist()
    # Restrict to expressed target genes
    tf_targets = tf_targets.map(lambda x: sorted(list(set(x) & set(gene))))
    # Restrict to TFs with at least min_targets genes
    targets =  tf_targets[tf_targets.map(lambda x: len(x) >= min_targets)]
    # print(targets)

    def fun_dca(v):
        global my_counter
        my_counter += 1
        print(f'{my_counter} / {len(targets)}')

        tmp = count[:,]
        tmp = ad.AnnData(tmp.X + 1)
        sc.pp.normalize_per_cell(tmp)
        size_factors = tmp.obs.n_counts/np.median(tmp.obs.n_counts)

        tmp = count[:,v]
        tmp = ad.AnnData(tmp.X + 1)
        tmp.obs["size_factors"]=size_factors

        dca(tmp, mode = 'latent', ae_type = ae_type, epochs=epochs,
                       early_stop=early_stop, hidden_size=hidden_size, verbose=verbose)
        return(tmp.obsm["X_dca"])
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
    if datatype == "continuous":
        rele_score = embed.map(fun_rfr)
    if datatype == "categorical":
        rele_score = embed.map(fun_rfc)
    return embed,rele_score

def rank_plot(result,save):
    score = pd.DataFrame(list(result[1].items()), columns=['Signature', 'Relevance Score'])
    score = score.sort_values('Relevance Score',ascending=False)
    new_df = score.head(5)
    new_df = new_df.append(score.tail(5))
    plt.figure(figsize=(10,10))
    ax = sns.barplot(x="Relevance Score", y="Signature", data=new_df, palette=sns.color_palette("Blues_r",10))
    plt.xticks(rotation='horizontal',fontsize=20)
    plt.yticks(fontsize=26)
    plt.xlabel('Relevance Score', fontsize=22)
    ax.grid(b=True, which='major', color='#d3d3d3', linewidth=0.5)
    fig = plt.gcf()
    plt.show()
    if save:
        fig.savefig(path + '/rank.svg', bbox_inches='tight')
    top_TF = score.head(n=5)['Signature'].tolist()
    print("Top_TF",":",top_TF)
    worse_TF = score.tail(n=5)['Signature'].tolist()
    print("Worse_TF",":",worse_TF)

def embedding_plot(result, tf_name, pheno, datatype, save):
    em = pd.DataFrame(result[0][tf_name],columns=['dca1','dca2'])
    plt.figure(figsize=(10, 8))
    if datatype == "continuous":
        plt.scatter(em.dca1, em.dca2,c = pheno)
        cbar = plt.colorbar()
        cbar.set_label("pseudotime", labelpad=+1)
    else:
        em['group'] = pheno.values
        em = em.sample(frac=1)
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x='dca1', y='dca2', hue="group", data=em, s=50)
    plt.title(tf_name, fontsize=30)
    plt.xlabel("dca1", fontsize=30)
    plt.ylabel("dca2", fontsize=30)
    fig = plt.gcf()
    plt.show()
    if save:
        fig.savefig(path + tf_name + '.svg', bbox_inches='tight')

def gene_plot(result, tf_name, gene, count, save):
    embedding = result[0][tf_name]
    plt.figure(figsize=(10,8))
    plt.title(gene,fontsize=30)
    plt.xlabel("dca1",fontsize=30)
    plt.ylabel("dca2",fontsize=30)
    sc.pp.scale(count)
    expr = count[:,gene].X
    expr[expr > 2] = 2
    plt.scatter(embedding[:, 0], embedding[:, 1],c = expr)
    cbar = plt.colorbar()
    cbar.set_label("expression", labelpad=+1)
    fig = plt.gcf()
    plt.show()
    if save:
        fig.savefig(path + gene + '.svg', bbox_inches='tight')
