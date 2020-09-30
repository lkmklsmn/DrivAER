import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numbers
import scanpy as sc
import rpy2.robjects as ro
import statsmodels.api as sm
from dca.api import dca

def rank_plot(result, save=False, path='./'):
    score = pd.DataFrame({'Signature': result[2], 'Relevance Score': result[1]})
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
        fig.savefig(path + 'rank.svg', bbox_inches='tight')
    top_TF = score.head(n=5)['Signature'].tolist()
    print("Top_TF",":",top_TF)
    worse_TF = score.tail(n=5)['Signature'].tolist()
    print("Worse_TF",":",worse_TF)

def embedding_plot(result, tf_name, pheno, save = False, path='./'):
    res = {result[2][i]: result[0][i] for i in range(len(result[2]))}
    em = pd.DataFrame(res[tf_name],columns=['dca1','dca2'])
    plt.figure(figsize=(10, 8))

    if isinstance(pheno[0], numbers.Number):
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

def gene_plot(result,tf_name,gene,count,pheno,save = False, path='./'):
    res = {result[2][i]: result[0][i] for i in range(len(result[2]))}
    em = pd.DataFrame(res[tf_name],columns=['dca1','dca2'])
    plt.figure(figsize=(10,8))
    if isinstance(pheno[0], numbers.Number):
        sc.pp.scale(count)
        expr = count[:,gene].X
        expr[expr > 2] = 2
        plt.scatter(em['dca1'], em['dca2'],c = expr)
        cbar = plt.colorbar()
        cbar.set_label("expression", labelpad=+1)
    else:
        em = em.set_index(count.obs.index)
        em['Value'] = pd.DataFrame(count[:,gene].X).values
        plt.figure(figsize=(10,8))
        cmap = sns.cubehelix_palette(dark=.2, light=0.9, as_cmap=True)
        sns.scatterplot(x='dca1', y='dca2', hue="Value",hue_norm=(-2,2),data=em,s=30,palette=cmap)
    plt.title(gene,fontsize=30)
    plt.xlabel("dca1",fontsize=30)
    plt.ylabel("dca2",fontsize=30)
    fig = plt.gcf()
    plt.show()
    if save:
        fig.savefig(path + gene + '.svg', bbox_inches='tight')
