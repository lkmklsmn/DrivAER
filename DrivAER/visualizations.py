import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numbers
import scanpy as sc
import rpy2.robjects as ro

def rank_plot(result, save = False):
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

def embedding_plot(result, tf_name, pheno, save = False):
    em = pd.DataFrame(result[0][tf_name],columns=['dca1','dca2'])
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

def gene_plot(result, tf_name, gene, count, save = False):
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

def heatmap():
    str="""
    plot(0)
    """
    ro.r(str)

def heatmap_pt(gene_set,tf_name,count, pseudo):
    targets = gene_set[tf_name]
    targets = list(set(targets) & set(count.var_names))
    c_o = count[:,targets]
    c = pd.DataFrame(c_o.X, columns=targets, index = count.obs.index.tolist())
    c['pseudotime'] = pseudo.values
    c = c.sort_values(by=['pseudotime'])
    # loess_fit
    def loess_fit(x):
        fit = sm.nonparametric.lowess(x.values, c['pseudotime'].values, frac=.6,it=0,is_sorted=True)[:,1]
        return (fit-fit.mean())/fit.std()
    fitted = c.apply(loess_fit, axis=0).drop('pseudotime',axis=1).T
    # time
    t = pseudo.sort_values(ascending=True)
    t = pd.DataFrame({'Pseudotime':t.values}).transpose()
    # TF dca coordinates
    coord = pd.DataFrame(result[0][tf_name],columns=['dca1','dca2'])
    coord['pseudotime'] = pseudo.values
    coord = coord.sort_values(by=['pseudotime']).transpose()
    coord = coord.drop('pseudotime',axis=0)
    # heatmap
    fig = plt.figure(figsize = (20,15))
    ax1 = plt.subplot2grid((20,20), (0,0), colspan=19, rowspan=1)
    ax1.set_title(tf_name)
    ax2 = plt.subplot2grid((20,20), (1,0), colspan=19, rowspan=2)
    ax3 = plt.subplot2grid((20,20), (3,0), colspan=19, rowspan=18)
    sns.heatmap(t,ax=ax1,xticklabels=False,cmap="rocket_r")
    sns.heatmap(coord,ax=ax2,xticklabels=False,cmap="rocket_r")
    ax1.set_yticklabels(labels=['pseudotime'],rotation=0)
    ax2.set_yticklabels(labels=['TF_dca'],rotation=0)
    sns.heatmap(fitted, ax=ax3, cmap="YlGnBu",xticklabels=False)
    
def heatmap_group(gene_set,tf_name,count,group,result):
    targets = gene_set[tf_name]
    targets = list(set(targets) & set(count.var_names))
    c_o = count[:,targets]
    c = pd.DataFrame(c_o.X, columns=targets, index = count.obs.index.tolist()).T
    fit = c.apply(lambda x:(x-x.mean())/x.std(),axis=1)              
    # time
    t = pd.DataFrame({'Cluster':pd.factorize(group)[0]}).transpose()
    # TF dca coordinates
    tmp = result[0][tf_name]
    coord = pd.DataFrame(tmp,columns=['dca1','dca2']).transpose()
    # heatmap
    fig = plt.figure(figsize = (20,15))
    ax1 = plt.subplot2grid((20,20), (0,0), colspan=19, rowspan=1)
    ax1.set_title(tf_name)
    ax2 = plt.subplot2grid((20,20), (1,0), colspan=19, rowspan=2)
    ax3 = plt.subplot2grid((20,20), (3,0), colspan=19, rowspan=18)
    sns.heatmap(t,ax=ax1,xticklabels=False,cmap="rocket_r")
    sns.heatmap(coord,ax=ax2,xticklabels=False,cmap="rocket_r")
    ax1.set_yticklabels(labels=['Cluster'],rotation=0)
    ax2.set_yticklabels(labels=['TF_dca'],rotation=0)
    sns.heatmap(fit, ax=ax3, cmap="YlGnBu",xticklabels=False,vmax=2)#,col_cluster=False)
