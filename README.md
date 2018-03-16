# t-SNE Heatmaps

## Introduction
When exploring a scRNA-seq dataset, it is common to visualize the cells with t-SNE and color them by expression of different genes that are markers for known cell types. The expression patterns of these marker genes suggest which clusters correspond to which known cell types. However, the expression of these marker genes can only be visualized one at a time, and only so many can fit on a page and compared side-by-side. We present t-SNE Heatmaps, which use 1D t-SNE to address this constraint, allowing for visualization of the expression patterns of hundreds of genes simulatneously.

## Example
In `1d_tsne_heatmaps_tutorial.R` we demonstrate t-SNE heatmaps using the sCNRA-seq data from Macosko et al. (2015) [[1]](http://www.cell.com/abstract/S0092-8674(15)00549-8).

Here is a screenshot of the result:
![alt text](http://gauss.math.yale.edu/~gcl22/t-sne-heatmaps/Macosko_etal "Example t-SNE heatmap")
**The interactive version is available [here](http://gauss.math.yale.edu/~gcl22/t-sne-heatmaps/Macosko_etal.html)**, and thanks to heatmaply, allows for more detailed exploration by zooming, panning, and more.

## How it works
Linderman and Steinerberger (2017) [[2]](https://arxiv.org/abs/1706.02582) prove that t-SNE will faithfully embed well clustered data, independent of the dimension. Real life data is not well clustered, but in our experience,  2D embeddings (in particular the clusters) tend to be generally consistent with 1D embedding.  However, 1D embeddings are much more compact, and allow us to build the heatmap visualization. To build a t-SNE heatmap, we do the following:

1. Compute 1D t-SNE of the cells. With FIt-SNE [[3]](https://arxiv.org/abs/1712.09005), this can be done for millions of cells in very short amount of time.
2. Discretize the 1D t-SNE embedding into 100 bins.
3. Represent each gene by the sum of its expression in the cells contained in each bin. Each gene is now a vector in R^100.
4. Visualize these vectors in heatmap format (i.e. each row is a gene and each column is a bin) using heatmaply.

The representation of the genes by their binned expression pattern in 1D t-SNE induces a distance function on the genes. In particular, we can perform hierarchical clustering on the genes using this new distance function, which often gives more meaningful results than Euclidean distance on the original data. This distance function also allows us to find genes whose expression pattern is associated with genes of interest. The user provides a list of genes of interest, and the algorithm then "enriches" this set with genes that have a similar expression pattern in the t-SNE.


## References
1. Macosko, Evan Z., et al. Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. Cell 161.5 (2015): 1202-1214.
2. Linderman, G. C., & Steinerberger, S. (2017). Clustering with t-SNE, Provably. arXiv:1706.02582
3. George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger, Yuval Kluger. (2017). Efficient Algorithms for t-distributed Stochastic Neighborhood Embedding. (2017) arXiv:1712.09005
