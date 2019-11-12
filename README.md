# DrivAER for manifold interpretation of scRNA-seq data
**DrivAER** is a method for identification of **Driv**ing transcriptional programs based on **A**uto**E**ncoder derived **R**elevance scores. 
It infers relevance scores for transcriptional programs with respect to specified outcomes of interest, which allows researchers to  identify when and where transcriptional programs are being up/down-regulated with high cellular resolution.

See our manuscript and [tutorial](https://github.com/lkmklsmn/TFscoring/blob/master/DrivAER_Tutorial.ipynb) for more details.

## Installation
### pip
	pip install -i https://test.pypi.org/simple/ DrivAER==0.0.1

## Input
1. raw count matrix
2. outcome of interest (pseudotemporal ordering/cell grouping etc)
3. gene set annotation

## Results
1. Relevance scores for each annotated transcriptional program
2. Data manifolds derived from each transcriptional program
3. Various visualizations (heatmap, DCA embedding, barplots)

## Usage

### Load annotation 
#### for gene set annotations in gmt format
| Set | Source | Target1 | Tatget2 | Target3|
| ---------- | ---------- |  :----:  |  :----:  |  :----:  | 
| set1 | source | gene1 | gene2 | gene3 |
| set2 | source | gene1 | gene2 | gene3 |
| set3 | source | gene1 | gene2 | gene3 |
	C3_mouse = get_anno(filename="C3.gmt",filetype="gmt",,conv_mouse=True)
#### for get set pairs in tsv format
| Set | Target | Type | Source|
| ---------- | ---------- |  :----:  |  :----:  | 
| set1 | gene1 | XX | XX |
| set1 | gene2 | XX | XX |
| set1 | gene3 | XX | XX |
| set2 | gene1 | XX | XX |
	trrust_human = get_anno(filename="trrust_human.tsv",filetype="tsv",conv_mouse=False)
### Calculate relevance scores
	import DrivAER as dv
	res = dv.calc_relevance(count = your_count, pheno = your_pt, tf_targets = C3_mouse, min_targets=5, datatype = "continuous")
### Generate visualizations
	dv.rank_plot(res, save)
	dv.embedding_plot(result, tf_name, pheno, datatype, save)
	dv.gene_plot(result, count, tf_name, gene, save)
	dv.heatmap(result, tf_name, save)
