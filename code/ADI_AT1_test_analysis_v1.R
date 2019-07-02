# Load R libs ####
library(plotrix)
library(mgcv)
library(slingshot)
library(SingleCellExperiment)
library(mclust)
library(RColorBrewer)
library(reshape2)
library(gam)
library(pheatmap)

load("../data/TFscoring_test.RData")

# Fit trajectory using slingshot ####
set.seed(123)
sim <- SingleCellExperiment(assays = List(counts = cm))
reducedDims(sim) <- SimpleList(DiffMap = data.matrix(dm))
cl1 <- Mclust(dm, G = 5)$classification
colData(sim)$GMM <- cl1
sce <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'DiffMap')

# Define pseudotime ####
t <- sce$slingPseudotime_1
t <- (t - min(t))/(max(t) - min(t))
if(cor(t, cm["Ager",]) < 0) t <- 1 - t

# Generate plots of embedding ####
par(mfrow = c(1, 3))
plot(reducedDims(sce)$DiffMap, pch=16, asp = 1, col = brewer.pal(9,"Set1")[cl1])
lines(SlingshotDataSet(sce), lwd=2, col='black')
legend('topleft', as.character(unique(cl1)), col = brewer.pal(9,"Set1")[unique(cl1)], pch = 16, bty = 'n')
plot(reducedDims(sce)$DiffMap, pch=16, asp = 1, col = color.scale(t, extremes = c('red', 'blue')))
boxplot(split(t, cl1)[names(sort(unlist(lapply(split(t, cl1), median))))], ylab = 'Pseudotime')

# Subset to genes expressed in at least 5 cells in more than 5 mice ####
cm <- assays(sce)$counts
md <- md[colnames(cm),]
genes <- rownames(cm)[which(apply(cm, 1, function(x) sum(x > 0)) > ncol(cm)/20)]
expr <- data.matrix(cm[genes,])

# Load TF annotation ####
tfs <- read.delim("../data/trrust_rawdata.mouse.tsv", header = F)

# Restrict to TFs present in the expression data with more than 5 target genes ####
ok <- intersect(tfs$V1, rownames(cm))
targets <- lapply(ok, function(x){
  intersect(tfs$V2[which(tfs$V1 == x)], rownames(cm))
})
names(targets) <- ok
targets <- targets[which(unlist(lapply(targets, length)) > 5)]

# Get TF scores ####
getTFscores <- function(expr = cm, targets = targets, pt = t){
  require(zoo)
  expr <- expr[,order(pt)]
  pt <- pt[order(pt)]
  output <- do.call(rbind, lapply(names(targets), function(nom){
    #print(nom)
    genes <- targets[[nom]]
    calcProb <- function(x){
      ok <- which(pt > min(x) & pt < max(x))
      
      global_prop <- apply(data.matrix(expr[genes,]), 1, function(y) sum(y > 0))
      global_prop <- data.frame(global_prop, ncol(expr))
      
      matr <- data.matrix(expr[genes, ok])
      matr <- matr[which(apply(matr, 1, function(y) sum(y > 0) > 0.1*ncol(matr))),]
      if(is.null(dim(matr))) return(NA)
      if(nrow(matr) < 3) return(NA)
      
      local_prop <- apply(matr, 1, function(y) sum(y > 0))
      local_prop <- data.frame(local_prop, ncol(matr))
      
      genes_ok <- intersect(rownames(local_prop), rownames(global_prop))
      
      if(length(genes_ok) < 3) return(NA)
      
      aframe <- data.frame(local_prop[genes_ok,], global_prop[genes_ok,])
      
      pvals <- unlist(lapply(1:nrow(aframe), function(y){
        prop.test(c(aframe[y, 1], aframe[y, 3]), c(aframe[y, 2], aframe[y, 4]), alternative = "greater")$p.value 
      }))
      -log10(fisher.method(t(data.matrix(pvals)))[,"p.value"])
    }
    rollapply(pt, width = 100, by = 25, FUN = calcProb, align = "left")
  }))
  rownames(output) <- names(targets)
  #output <- t(apply(output, 1, function(x) 1 - (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))))
  output[which(is.na(output))] <- 0
  output
}
res <- getTFscores(expr = cm, targets = targets, pt = t)
empty <- which(apply(res, 1, var) == 0)
if(length(empty) > 0) res <- res[-empty,]
res[which(res > 10)] <- 10

# Run loess regression ####
loess_res <- apply(res, 1, function(x){
  aframe <- data.frame(activity = x, pt = 1:ncol(res))
  aframe$activity[which(is.na(aframe$activity))] <- 0
  gam(activity ~ lo(pt), data = aframe)
})
loess_pvals <- unlist(lapply(loess_res, function(x) summary(x)[4][[1]][1,5]))

# Generate heatmap of significant activity scores ####
sig <- which(p.adjust(loess_pvals, method = "BH") < 0.25)
loess_fitted <- do.call(rbind, lapply(loess_res[sig], function(x) predict(x, data.frame(pt = 1:ncol(res)))))
tmp <- loess_fitted
peaks <- apply(tmp, 1, function(x) which(x == max(x)))
rowOrd <- order(peaks)
anno_col <- data.frame(pt = seq(min(t), max(t), length = 50))
rownames(anno_col) <- colnames(tmp)
pheatmap(tmp[rowOrd,], scale = "row", cluster_cols = F, cluster_rows = F, annotation_col = anno_col, show_colnames = F)

# Define additional plotting functions ####
plotScore <- function(nom){
  scale <- function(x) (x - mean(x))/sd(x)
  #zscore <- t(apply(log(cm[targets[[nom]],] + 1), 1, scale))
  zscore <- apply(cm[targets[[nom]],], 2, function(x) mean(x > 0))
  #zscore <- apply(zscore, 2, function(x) mean(x, na.rm = T))
  #zscore <- scale(zscore)
  #zscore[which(zscore > 2)] <- 2
  zscore <- scale(zscore)
  farben <- color.scale(zscore, extremes = c("blue", "grey", "red"), alpha = 0.8)
  #plot(dm, col = farben, pch = 16, main = nom)
  aframe <- data.frame(t, zscore)
  ggplot(aframe, aes(y = zscore, x = t)) +
    geom_smooth(method = "loess") +
    ylab("TF activity score") + xlab("Pseudotime") + ggtitle(nom)
}
plotGenes <- function(nom){
  genes <- targets[[nom]]
  expr <- log(data.matrix(cm[genes,]) + 1)
  expr <- expr[which(apply(expr, 1, var) > 0),]
  expr <- t(apply(expr, 1, function(x) (x - mean(x))/sd(x)))
  aframe <- data.frame(pt = t, t(expr))
  melted <- melt(aframe, id.vars = "pt", measure.vars = setdiff(colnames(aframe), "pt"))
  if(length(genes) > 10){
    ggplot(melted, aes(y = value, x = pt, color = variable)) +
      geom_smooth(method = "loess", span = 1, show.legend = F) +
      ylab("Relative expression") + xlab("Pseudotime") + ggtitle(nom)  
  }
  else{
    ggplot(melted, aes(y = value, x = pt, color = variable)) +
      geom_smooth(method = "loess", span = 1) +
      ylab("Relative expression") + xlab("Pseudotime") + ggtitle(nom)
  }
}

# Plot some examples ####
plotScore("Trp53")
plotGenes("Sox17")
