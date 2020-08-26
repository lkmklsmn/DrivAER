# Simulating path
library(splatter)
library(scater)
params.groups <- newSplatParams(batchCells = 2000, nGenes = 500)
sim <- splatSimulatePaths(params.groups,
                          group.prob = c(1),
                          de.prob = 0.2, de.facLoc = 0.2,
                          path.from = c(0),
                          dropout.type = 'experiment',
                          dropout.shape=-1,dropout.mid=1, 
                          verbose = FALSE)
sim <- logNormCounts(sim)
sim <- runPCA(sim)
sim <- runUMAP(sim,dimred="PCA",n_dimred=10)
plotReducedDim(sim, dimred = "UMAP", colour_by = "Step")+ ggtitle("UMAP with dropout, de.prob=0.2, 5 PCs")

