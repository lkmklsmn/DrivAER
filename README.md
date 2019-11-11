# DrivAER for Manifold interpretation of scRNA-seq data
DrivAER is a method for identification of **Driv**ing transcription programs based on **A**uto**E**ncoder derived **R**elevance scores. 
It infers relevance scores for transcription programs with respect to specified variables of interest, which allows researchers to prioritize transcription programs and identify when and where transcription programs are being up/down-regulated with high cellular resolution.

See our manuscript and [tutorial](https://github.com/lkmklsmn/TFscoring/blob/master/DrivAER_Tutorial.ipynb) for more details.

## Installation
### pip
	pip install -i https://test.pypi.org/simple/ DrivAER==0.0.1

## Input
1. raw count matrix
2. outcome of interest (pseudotemporal ordering/cell grouping)
3. (optional) user-defiend transcriptional gene set

## Results
1. Ranking plot of relevance scores for annotated transcription programs
2. Data manifold derived from a specific transcriptional program
3. Heatmap of a specific transcriptional program

## Usage
### Load raw count matrix and phenotype
	your_count
	your_pt
### Load target gene set 
	C3_mouse =get_anno(filename="C3.gmt",transfer=True)
### Run enrich_test
	import DrivAER as dv
	res = dv.calc_relevance(count = your_count, pheno = your_pt, tf_targets = C3_mouse, min_targets=5, datatype = "continuous")
### Get output plot
	dv.rank_plot(res, save)
	dv.embedding_plot(result, tf_name, pheno, datatype, save)
	dv.gene_plot(result, count, tf_name, gene, save)
	dv.heatmap(result, tf_name, save)
