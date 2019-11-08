# DrivAER for Manifold interpretation
DrivAER (identification of **Driv**ing transcription programs based on **A**uto**E**ncoder derived **R**elevance scores) infers relevance scores for transcription programs with respect to specified variables of interest, which allows researchers to prioritize transcription programs and identify when and where transcription programs are being up/down-regulated with high cellular resolution.

See our manuscript and tutorial for more details.

## Installation
### pip
	pip install -i https://test.pypi.org/simple/ DrivAER==0.0.1

## Input
1. raw count matrix
2. variables of interest (cluster/pseudotemporal trajectory)
3. (optional) user-defiend target gene set

## Results
1. Ranking plot of relevance scores for annotated transcription programs
2. Data manifold derived from a specific transcription program
3. Heatmap of a specific transcriptional program

## Usage
### raw count matrix
	count = sc.read("file.txt", sep="\t", first_column_names=True,cache=True)
### pseudotime
	pt = pd.read_csv("file.txt")
### user-defined target gene set 
	C3_mouse =get_anno(filename="C3.gmt",transfer=True)
### Run enrich_test
	res = enrich_test(count = count, pheno = pt, tf_targets = C3_mouse, min_targets=5, datatype = "continuous")
### Get output plot
	rank_plot(res,save=True)
	enrich_plot_tf(result, tf_name, pheno, datatype, save)
	enrich_plot_marker(result, tf_name, gene, count, save)
	heatmap(result,tf_name)
