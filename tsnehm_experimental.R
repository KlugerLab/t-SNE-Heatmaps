#' Compute t-SNE Heatmap
#'
#' Experimental version of t-SNE heatmaps codes
#'
#' @param expression_matrix Full expression matrix genes (rows) vs. cells (columns). Rownames should be gene names
#' @param goi goi Genes of Interest
#' @param tsne_embedding 1D t-SNE embedding
#' @param cell_labels Labels for each cell. Typically these are the cluster assignments for each cell, as obtained by 
#'                    dbscan on the 2D t-SNE. This way, the columns of the heatmap will have a color assigned to each of them,
#'                    and they can be mapped to a corresponding location on the 2D t-SNE
#' @param enrich For every gene in goi, find enrich number of genes that are close to that gene in the distance induced by the 1D t-SNE
#' @param breaks Number of bins
library(heatmaply)
require(gplots)
library(pdist)
library(plyr)
library(RColorBrewer)
library(parallel)
library(data.table)

genes_tsne_vecs <- function(expression_matrix, tsne_embedding,  breaks=100 ){
  
  nonempty_cells <- colSums(expression_matrix) > 0;
  expression_matrix <- expression_matrix[,nonempty_cells]
  tsne_embedding <- tsne_embedding[nonempty_cells]
  
  nonempty_genes <- rowSums(expression_matrix) > 0;
  expression_matrix <- expression_matrix[nonempty_genes,]
  
  goidf <- data.frame(x=tsne_embedding, t(expression_matrix))
  group <- split (goidf, cut(goidf$x,breaks = breaks))
  apply_step <- mclapply(group,function(x) as.data.frame(t(as.data.frame(colSums((x[,-1])))),mc.cores = 10))
  bin_counts <- as.data.frame(data.table::rbindlist(apply_step))
  
  outmat <- sweep(bin_counts,2,colSums(bin_counts), '/')
  colnames(outmat) <- rownames(expression_matrix)
  outmat
}
tsnehm <- function(expression_matrix, goi, tsne_embedding, cell_labels, enrich=0, breaks=100, slope=50, intercept=0.05){
  
  nonempty_cells <- colSums(expression_matrix) > 0;
  expression_matrix <- expression_matrix[,nonempty_cells]
  tsne_embedding <- tsne_embedding[nonempty_cells]
  cell_labels <- cell_labels[nonempty_cells]
  
  nonempty_genes <- rowSums(expression_matrix) > 0;
  expression_matrix <- expression_matrix[nonempty_genes,]
  
  notinrows <- !(goi%in% rownames(expression_matrix));
  if (sum(notinrows) > 0){
    print(sprintf("The following rows are not in the matrix: %s", paste(goi[notinrows])))
  }
  goi <- goi[!notinrows]
  if (enrich ==0 ) {
    goidf <- data.frame(x=tsne_embedding, t(expression_matrix[goi,]))
  }else{
    goidf <- data.frame(x=tsne_embedding, t(expression_matrix))
  }
  group <- split (goidf, cut(goidf$x,breaks = breaks))
  apply_step <- mclapply(group,function(x) as.data.frame(t(as.data.frame(colSums((x[,-1])))),mc.cores = 10))
  bin_counts <- as.data.frame(data.table::rbindlist(apply_step))
  bin_counts_s  <- sweep(bin_counts,2,colSums(bin_counts), '/')
  if (enrich >0 ) {
    enriched_genes_list <- list();
    print("Now enriching")
    pdisttest <- pdist(t(as.matrix(bin_counts_s)),indices.A = goi, indices.B=1:ncol(bin_counts_s))
    sortedpdist <- t(apply(as.matrix(pdisttest), 1, order,decreasing=FALSE))
    enriched_genes <- unique(as.vector(t(sortedpdist[,1:(enrich +1)])))
    for (goii in 1:length(goi)){
      enriched_genes_list[[goi[goii]]] <-rownames(expression_matrix)[sortedpdist[goii,2:enrich+1] ]
    }
    bin_counts_s <- bin_counts_s[,enriched_genes]
    
    
  }
  
  #assign a label to each column based on which of the cell_labels is the most common
  dbscan_tsne <- data.frame(x=tsne_embedding, y=cell_labels)
  dbscan_group <- split (dbscan_tsne, cut(dbscan_tsne$x,breaks = breaks))
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  group_labels <-unlist(lapply(dbscan_group, function(x) Mode(x[,2])))
  rownames(bin_counts_s) <- sprintf("bin:%s, label: %d", rownames(bin_counts_s), group_labels)
  group_labels[is.na(group_labels)] <- NA
  if (enrich>0){
    gene_colors <- rep(0, length(enriched_genes));
  }else{
    gene_colors <-c();
  }
  #gene_colors[seq(1,length(enriched_genes), by=(enrich+1))] <- 1
  gene_colors[colnames(bin_counts_s) %in% goi] <- 1
  
 my_palette <- colorRampPalette(c("white", "red"))(n = 1000)
 row_color_palette <- colorRampPalette(c("white", "blue"))
 col_color_palette <- colorRampPalette(brewer.pal(n = max(group_labels+1,na.rm=TRUE),"Spectral"))
 toplot <- t(as.matrix(bin_counts_s,2))
 
 #toplot <- t(as.matrix(bin_counts_s>0.05,2) + 0)
 #
 #fofo <- seq(0,0.3,by = 0.001)
 #plot(fofo,1/(1+exp(50*0.01+-50*fofo)))
 #fifi <- apply(toplot,FUN = quantile,MARGIN = 1,0.99) 
 #toplot15 <- sweep(toplot,1,fifi/50 ,  "+")
 #toplot2 <- 1/(1+exp(50*fifi+-50*toplot15))
 #toplot15 <- sweep(toplot,2,fifi ,  "+")
  
  
 #Best so far
 #toplot2 <- 1/(1+exp(50*0.05-50*toplot))
 toplot2 <- 1/(1+exp(slope*intercept-slope*toplot))
  return_list <- list()
  return_list$heatmap <- heatmaply(toplot2,col_side_colors = group_labels, row_side_colors=gene_colors, 
                                   row_side_palette=row_color_palette, showticklabels = c(FALSE,TRUE),  
                                   col=my_palette,  dendrogram='row', titleX =FALSE, RowV=FALSE,
                                   show_legend=FALSE,hide_colorbar = TRUE,fontsize_row = 5,margins = c(70,50,NA,0)
  )
  if (enrich>0){
    
  return_list$enriched_genes <- enriched_genes_list
  }
  return_list
}