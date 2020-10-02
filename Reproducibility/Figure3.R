# Load R libs ####
library(ggplot2)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())

# Load accuracies ####
accs <- read.csv('data/accuracies.csv', row.names = 1)
tmp <- melt(accs, measure.vars = c('drivaer', 'pca', 'umap', 'tsne'))

# Plot Fig 3d #### 
ggplot(tmp, aes(x = fraction, y = value, group = variable, color = variable)) +
  geom_smooth() + ylab('Relevance score')

# Plot Fig 3g ####
tmp <- accs[which(accs$fraction %in% c(0,0.2)), ]
t.test(split(tmp$drivaer, tmp$fraction)[[1]], split(tmp$drivaer, tmp$fraction)[[2]],
       alternative = 'less')
ggplot(tmp, aes(x = fraction, y = drivaer, group = fraction, color = as.factor(fraction))) +
  geom_boxplot() + geom_point(color = 'black') + ylab('Relevance score')

# Load data from different bottleneck configs ####
accs <- read.csv('data/hidden_size_diff_accuracies.csv', row.names = 1)
tmp <- melt(accs, measure.vars = setdiff(colnames(accs), 'fraction'))

# Plot Fig 3f ####
ggplot(tmp, aes(x = fraction, y = value, group = variable, color = variable)) +
  geom_smooth() + ylab('Relevance score')

# Load accuracies from PAGODA and VISION ####
vision <- read.csv('data/VISION_autocorrelation.csv', row.names = 1)
pagoda <- read.csv('data/PAGADA_accuracies.csv', row.names = 1)
accs <- read.csv('data/accuracies.csv', row.names = 1)
tmp <- cbind(vision, pagoda = pagoda$adj.z, drivaer = accs$drivaer)
tmp <- melt(tmp, measure.vars = setdiff(colnames(tmp), 'fraction'))

# Plot Fig 3h,i,j ####
p1 <- ggplot(tmp[which(tmp$variable == 'pagoda'),], aes(x = fraction, y = value, group = variable, color = variable)) +
  geom_smooth() + theme(legend.position = "none") +
  scale_x_continuous(name="Fraction DE", breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(name='Adjusted z-score', limits=c(-3, 3))

p2 <- ggplot(tmp[which(tmp$variable %in% c('score_directed', 'score_undirected')),], aes(x = fraction, y = value, group = variable, color = variable)) +
  geom_smooth() +  theme(legend.position = "none") +
  scale_x_continuous(name="Fraction DE", breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(name="Autocorrelation", limits=c(0, 1))

p3 <- ggplot(tmp[which(tmp$variable == 'drivaer'),], aes(x = fraction, y = value, group = variable, color = variable)) +
  geom_smooth() + theme(legend.position = "none") +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_x_continuous(name="Fraction DE", breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(name="Relevance score", limits=c(0, 1))

gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
