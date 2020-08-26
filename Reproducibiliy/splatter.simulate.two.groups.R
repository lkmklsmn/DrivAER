# Load R libs ####
library(splatter)
library(scater)
library(Matrix)
library(Seurat)

# Simulate scRNAseq data ###
sim <- splatSimulate(group.prob=c(0.5,0.5), 
                     batchCells=2000,
                     nGenes=500,
                     de.prob= 0.1,
                     de.facLoc = 0.01,
                     dropout.type = 'experiment',
                     dropout.shape=-1,dropout.mid=5, 
                     method="groups")

# Plot PCAs ####
gene_info <- rowData(sim)
de <- rownames(gene_info)[which(gene_info$DEFacGroup1 != gene_info$DEFacGroup2)]
nonde <- setdiff(rownames(sim), de)

sim <- logNormCounts(sim)
sim <- runPCA(sim, ncomponents = 10)

p_all <- plotPCA(sim, colour_by = "Group") + ggtitle("All")
p_de <- plotPCA(runPCA(sim[de,], ncomponents = 10), colour_by = "Group") + ggtitle("DE only")
p_nonde <- plotPCA(runPCA(sim[nonde,], ncomponents = 10), colour_by = "Group") + ggtitle("non-DE only")
gridExtra::grid.arrange(p_all, p_de, p_nonde, ncol = 3)

# Save files ####
write.csv(colData(sim), file = "two_group_simulation_cell_info.csv")
write.csv(assays(sim)[['counts']], file = "two_group_simulation_counts.csv")
write.csv(gene_info, file = "two_group_simulation_gene_info.csv")