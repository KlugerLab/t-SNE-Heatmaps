# Download the data -----------------------------------------------------------
library(rsvd)
library(dbscan)
library(ggplot2)

matrix_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63472&format=file&file=GSE63472%5FP14Retina%5Fmerged%5Fdigital%5Fexpression%2Etxt%2Egz';
download.file(matrix_url,'MacoskoMatrix.txt.gz')
system('gunzip MacoskoMatrix.txt.gz')

A <- read.table("/data/Linderman/RWorking/MacoskoMatrix.txt",sep = "\t",header = TRUE)
rownames(A) <- A[,1]
A <- A[,-1]


# Preprocessing -----------------------------------------------------------
# Having loaded the data, let's do some (very) basic preprocessing. Divide by the
# number of UMIs per each cell (column), then multiply by 10E3.  In other words,
# for a given cell (column), each row is the % of the total UMIs contributed by
# these gene.   Then log normalize

totalUMIPerCell <- colSums(A);

A_norm <- sweep(A, 2, totalUMIPerCell, '/');
A_norm <- A_norm * 10E3
A_norm <- log(A_norm +1);

# PCA -----------------------------------------------------------
# Now let's do PCA. We will use the randomized SVD package rsvd. 
A_rm <- rowMeans(A_norm)
A_norm_c <- sweep(A_norm, 1,A_rm) 

fastDecomp <- rsvd(A_norm_c,k=50)
scores <- fastDecomp$v %*% diag(fastDecomp$d);
plot(scores[,1], scores[,2], pch='.')

# t-SNE -----------------------------------------------------------

# Obtain FIt-SNE from https://github.com/KlugerLab/FIt-SNE
setwd("/data/Linderman/FIt-SNE"); source('fast_tsne.R')
system.time(tsne.vptrees.1d <- fftRtsne(as.matrix(scores[,1:50]),  dims=1,ann_not_vptree = FALSE,
                                        data_path = "temp/data.dat", result_path = "temp/result.dat",rand_seed = 3))

#We don't need it for the heatmaps, but since we are used to looking at 2D t-SNE, let's compute that too 
system.time(tsne.vptrees.2d <- fftRtsne(as.matrix(scores[,1:50]),  dims=2,ann_not_vptree = FALSE,
                                        data_path = "temp/data.dat", result_path = "temp/result.dat",rand_seed = 3))

clustering <- dbscan(tsne.vptrees.2d,eps=1, minPts = 40)
ggplot(data.frame(V1=tsne.vptrees.2d[,1], V2=tsne.vptrees.2d[,2], color=as.factor(clustering$cluster)), aes(x=V2,y=V1, color=color) ) +  geom_point(size=0.05)  +
  guides(colour = guide_legend(override.aes = list(size=5), nrow=3)) + theme(legend.position="bottom",title = element_blank())

# Take a look at the 1D t-SNE, before we feed it into the heatmaps
# In order to see it better, let's make some noise for the Y axis
noise_y <- runif(length(tsne.vptrees.1d))
# Color it with a marker gene for a cluster from the paper to make sure it localizes
ggplot(data.frame(V1=noise_y, V2=tsne.vptrees.1d, 
                            col=as.factor(A_norm["LHX1",]>0)), aes(x=V2,y=V1, color=col) ) + geom_point(size=0.05)

# Make a t-SNE heatmap! -----------------------------------------------------------
source('/data/Linderman/tsne_applications/tsnehm_experimental.R')

# Here's a list of "genes of interest" that we think might mark a specific
# cluster. These are typically chosen by prior knowledge. Here, we pick genes
# that were used in Figure 5 of Macosko et al or Figure 1 of the related paper
# studying Bipolar cells (Shekhar et al. 2016).
querygoi <- c("LHX1", "SLC17A6",  "PAX6", "GAD1", "SLC6A9", "OPN1MW","VSX1",
         "RLBP1", "GFAP", "PECAM1", "KCNJ8", "CX3CR1",
         "VSX2", "OTX2", "SCGN", "ISL1", "GRM6", "APOE",  "RHO", "ARR3",
         "TACR3", "SYT2", "NETO1", "IRX6", "PRKAR2B", "GRIK1", "KCNG4", "CABP5", "PRKCA")

# This might take a little while. The reason is that we have so many genes. If 
# you are only interested in, say, the 1000 most variable genes, restrict the
# first argument to only those genes, and this step can be faster
hmtsneout <- tsnehm(A_norm, querygoi, tsne.vptrees.1d,clustering$cluster,enrich=5)
hmtsneout$heatmap
