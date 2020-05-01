library(splatter)
library(scater)
library(DEsingle)
load("lupus.sce.RData")

# Estimated parameters of stimulated group ####
stim <- rownames(colData(lupus.sce))[which(colData(lupus.sce)$label=="stim")]
subset <- lupus.sce[,stim]
param.stim <- splatEstimate(as.matrix(assay(subset)))
param.stim <- setParam(param.stim, "dropout.type", "experiment")
lupus.sim.stimu <- splatSimulate(param.stim,nGenes = 2000,batchCells=1000)

# Estimated parameters of control group ####
ctrl <- rownames(colData(lupus.sce))[which(colData(lupus.sce)$label=="ctrl")]
subset <- lupus.sce[,ctrl]
param.ctrl <- splatEstimate(as.matrix(assay(subset)))
param.ctrl <- setParam(param.ctrl, "dropout.type", "experiment")
lupus.sim.ctrl <- splatSimulate(param.ctrl,nGenes = 2000,batchCells=1000)

# Combined ####
lupus_sim <- cbind(counts(lupus.sim.stimu),counts(lupus.sim.ctrl))
group <- factor(c(rep(1,1000), rep(2,1000)))

# Find DE genes ####
results <- DEsingle(counts = lupus_sim, group = group)
results.classified <- DEtype(results = results, threshold = 0.05)
de.genes <- rownames(results.classified)[which(results.classified$Type=="DEg")]

# Change into SingleCellExperiment ####
lupus_sim <- SingleCellExperiment(lupus_sim)
colnames(lupus_sim) <- paste0("Cell",seq(1:2000))
rowData(lupus_sim)$label <- ifelse(rownames(rowData(sce)) %in% de.genes,"DE","non_DE")
save(param.stim,param.ctrl,lupus_sim,file="kang_simu_parameters.RData")
