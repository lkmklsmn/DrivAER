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
1. Ranking plot of relevance scores for annotated transcription programs
2. Data manifold derived from a specific transcriptional program
3. Heatmap of a specific transcriptional program

## Usage

### Load annotation 
#### for annotations in gmt format
	C3_mouse = get_anno(filename="C3.gmt",filetype="gmt",,conv_mouse=True)
#### for TF-target pairs in tsv format
	trrust_human = get_anno(filename="trrust_human.tsv",filetype="tsv",conv_mouse=False)
### Calculate relevance scores
	import DrivAER as dv
	res = dv.calc_relevance(count = your_count, pheno = your_pt, tf_targets = C3_mouse, min_targets=5, datatype = "continuous")
### Generate visualizations
	dv.rank_plot(res, save)
	dv.embedding_plot(result, tf_name, pheno, datatype, save)
	dv.gene_plot(result, count, tf_name, gene, save)
	dv.heatmap(result, tf_name, save)
