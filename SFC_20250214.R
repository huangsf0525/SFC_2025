
####====================  Program Description  ====================####

## The SFC (Segmentation Fusion Cluster) algorithm segments the features of the dataset, 
## constructs a similarity matrix, applies a fusion method, and finally calculates the 
## clustering index Ψ to evaluate the clustering performance.

####====================  Parameter Description  ====================####
## Input
## data   : Dataset, where the i-th row represents the i-th observation 
##          and the j-th column represents the j-th variable 
##          for i=1,..,n and j=1,...,p.
## m      : Number of segments for m > 1
## cluster: Number of clusters
## K      : The range of neighborhood in the fusion step. Default value is
##          0.6 X p/cluster 
##
## Output
## Ψ      : The proposed evaluation metrics of SFC_E1, SFC_1, and SFC_m 
##          with spectral, k-means, and hierarchical clustering methods


SFC <- function(data = NULL, m = NULL, cluster = NULL, K = floor(0.6*ncol(selected_data)/c)) {

  ####====================  function  ====================####
  library(SNFtool)
  library(philentropy)
  library(caret)
  library(ggplot2)
  library(cluster)
  library(clue)
  
  ### The "affinityMatrix" function constructs similarity networks
  affinityMatrix <- function(Diff,K=20,sigma=0.5) {
    N = nrow(Diff)
    
    Diff = (Diff + t(Diff)) / 2
    diag(Diff) = 0;
    sortedColumns = as.matrix(t(apply(Diff,2,sort)))
    finiteMean <- function(x) { mean(x[is.finite(x)]) }
    means = apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
    
    avg <- function(x,y) ((x+y)/2)
    Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
    Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
    densities = dnorm(Diff,0,sigma*Sig,log = FALSE)
    
    changedensities=densities*sqrt(2*pi)*Sig*sigma
    W = (changedensities + t(changedensities)) / 2
    
    return(W)
  }
  
  ### The "distrue" function computes the Ψ value
  distrue = function(mat_normalized,cluster){
    
    group            = cluster
    ordered_indices  = order(group)
    reordered_matrix = mat_normalized[ordered_indices, ordered_indices]
    
    # Within clusters
    C     = list()
    mean1 = c()
    l     = c()
    for (i in 1:length(unique(group))) {
      C[[i]]       = as.matrix(log(reordered_matrix[(length(which(group<=i-1))+1):length(which(group<=i)),(length(which(group<=i-1))+1):length(which(group<=i))]))
      diag(C[[i]]) = 0
      mean1[i]     = sum(C[[i]])
      l[i]         = length(C[[i]])
    }
    
    # Between clusters
    compute    = reordered_matrix
    logcompute = log(compute)
    for (i in 1:length(unique(group))) {
      logcompute[(length(which(group<=i-1))+1):length(which(group<=i)),(length(which(group<=i-1))+1):length(which(group<=i))] = 0  
    }
    logical_index = row(logcompute) < col(logcompute)
    result        = logcompute[logical_index]
    sum2          = sum(result)
    logical_index = row(logcompute) > col(logcompute)
    result        = logcompute[logical_index]
    sum3          = sum(result)
    
    d             = ((sum2+sum3)/(length(logcompute)-sum(l))) / (sum(mean1)/(sum(l)-ncol(mat_normalized)))
    return(d)}
  
  ### The "odistrue" function computes the Ψ value
  odistrue = function(data, cluster, K = floor(0.6*nrow(data)/length(unique(cluster)))){
    alpha                = 0.5;  	
    T                    = 20; 
    data                 = standardNormalization(data)
    Dist1                = (dist2(as.matrix(data),as.matrix(data)))^(1/2);
    w1                   = affinityMatrix(Dist1, K, alpha)
    w                    = w1
    diag(w1)             = 0
    row_sums             = rowSums(w1)
    mat_normalized       = (w1/ (2*row_sums))
    diag(mat_normalized) = 0.5
    group                = cluster
    ordered_indices      = order(group)
    reordered_matrix     = mat_normalized[ordered_indices, ordered_indices]
    
    # Within clusters
    C     = list()
    mean1 = c()
    l     = c()
    for (i in 1:length(unique(group))) {
      C[[i]]       = as.matrix(log(reordered_matrix[(length(which(group<=i-1))+1):length(which(group<=i)),(length(which(group<=i-1))+1):length(which(group<=i))]))
      diag(C[[i]]) = 0
      mean1[i]     = sum(C[[i]])
      l[i]         = length(C[[i]])
    }
    
    # Between clusters
    compute    = reordered_matrix
    logcompute = log(compute)
    for (i in 1:length(unique(group))) {
      logcompute[(length(which(group<=i-1))+1):length(which(group<=i)),(length(which(group<=i-1))+1):length(which(group<=i))] = 0  
    }
    logical_index = row(logcompute) < col(logcompute)
    result        = logcompute[logical_index]
    sum2          = sum(result)
    logical_index = row(logcompute) > col(logcompute)
    result        = logcompute[logical_index]
    sum3          = sum(result)
    
    d             = ((sum2+sum3)/(length(logcompute)-sum(l))) / (sum(mean1)/(sum(l)-ncol(mat_normalized)))
    return(list(d,w=w,cw=mat_normalized))}
  
## Main program ------------------
  
if (!is.matrix(data) || is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    cat("Error: 'data' must be a non-empty matrix.\n")
    return(NULL)
}
    
if (!is.numeric(cluster) || cluster <= 0 || cluster != as.integer(cluster) || is.na(cluster)) {
    cat("Error: 'cluster' must be a positive integer.\n")
    return(NULL)
}
    
selected_data <- data
c <- cluster
p <- dim(selected_data)[1]
max_p <- floor(p/2)
    
if (!is.numeric(m) || m <= 0 || m != as.integer(m) || is.na(m)) {
    cat("Error: 'm' must be a positive integer.\n")
    return(NULL)
} else if (m < 1 || m > max_p) {
    cat("Error: 'm' must be between 1 and", max_p, "\n")
    return(NULL)
}
    
all_results_acc <- list()
all_results_dis <- list()
    
# Method: SFC_E1
{
set.seed(123)
data = t(selected_data)
Test1_euclidean = suppressMessages(distance(data,method = "euclidean", p = NULL, test.na = TRUE, unit = "log")) 
diag(Test1_euclidean) = 0
row_sums <- rowSums(Test1_euclidean)
mat_normalized1 <- (Test1_euclidean/ (2*row_sums))
diag(mat_normalized1) = 0.5
        
#spectral clustering
expect_label_sp = spectralClustering(mat_normalized1,c) 
euclidean_dis_sp = distrue(mat_normalized1,expect_label_sp)
        
#kmeans
Test1_euclidean_kmeans = kmeans(Test1_euclidean,centers = c,nstart = 20) 
expect_label_kmeans = Test1_euclidean_kmeans$cluster 
euclidean_dis_kmeans = odistrue(data,expect_label_kmeans)[[1]]
        
#hierarchical clustering
Test1_euclidean_hc = hclust(dist(data, method = 'euclidean'),method = "ward.D2")
expect_label_hc = cutree(Test1_euclidean_hc,k = c) 
euclidean_dis_hc = odistrue(data,expect_label_hc)[[1]]
        
euclidean_result2 = data.frame(spec = euclidean_dis_sp, kmeans = euclidean_dis_kmeans, hc = euclidean_dis_hc)
rownames(euclidean_result2) <- "Ψ"
all_results_dis[["SFC_E1"]] <- euclidean_result2
}
    
# Method: SFC_1
{
data <- t(selected_data)  
data <- standardNormalization(data)                            
Dist <- sqrt(dist2(data, data))  
O <- affinityMatrix(Dist, K, 0.5) 
        
w1 <- O
diag(w1) <- 0  
row_sums <- rowSums(w1)  
mat_normalized2 <- (w1 / (2 * row_sums))  
diag(mat_normalized2) <- 0.5
        
# Spectral clustering
Olabel_sp <- spectralClustering(mat_normalized2, c)
SNF_affinity_dis_sp <- odistrue(data, Olabel_sp)[[1]]
        
# K-Means
kmeans_O <- kmeans(O, centers = c, nstart = 20)
expect_label_kmeans <- kmeans_O$cluster
SNF_affinity_dis_kmeans <- odistrue(data, expect_label_kmeans)[[1]]
        
# Hierarchical Clustering
hc_O <- hclust(dist(O, method = 'euclidean'), method = 'ward.D2')
expect_label_hc <- cutree(hc_O, c)
SNF_affinity_dis_hc <- odistrue(data, expect_label_hc)[[1]]
        
SNF_affinity_result2 <- data.frame(
        spec = SNF_affinity_dis_sp,
        kmeans = SNF_affinity_dis_kmeans,
        hc = SNF_affinity_dis_hc)
rownames(SNF_affinity_result2) <- "Ψ"

all_results_dis[["SFC_1"]] <- SNF_affinity_result2
}
    
# Method: SFC_m
if (m >1) {
   set.seed(111225008)
   random_nums <- sample(1:p, p)
   group_sizes <- rep(floor(p / m), m)
   if (p %% m > 0) {
      group_sizes[1:(p %% m)] <- group_sizes[1:(p %% m)] + 1
   }
   groups <- split(random_nums, rep(1:m, each = floor(p / m), length.out = p))
        
   affinity_matrices <- list()
   for (i in seq_along(groups)) {
       group_data <- selected_data[groups[[i]], ]
       group_data <- t(group_data)
       group_data <- standardNormalization(group_data)
       group_dist <- sqrt(dist2(group_data, group_data))
       affinity_matrices[[i]] <- affinityMatrix(group_dist, K, 0.5)
   }
        
   w <- SNF(affinity_matrices, K, 20)
        
   # Spectral Clustering
   expect_label_sp <- spectralClustering(w, c)
   split_fusion_dis_sp <- distrue(w, expect_label_sp)
       
   # K-Means
   expect_label_kmeans <- kmeans(w, centers = c, nstart = 20)$cluster
   split_fusion_dis_kmeans <- distrue(w, expect_label_kmeans)
        
   # Hierarchical Clustering
   hc_split_w <- hclust(dist(w, method = 'euclidean'), method = 'ward.D2')
   expect_label_hc <- cutree(hc_split_w, c)
   split_fusion_dis_hc <- distrue(w, expect_label_hc)
        
   split_fusion_result2 <- data.frame(spec = split_fusion_dis_sp, kmeans = split_fusion_dis_kmeans, hc = split_fusion_dis_hc)
   rownames(split_fusion_result2) <- "Ψ"
        
   all_results_dis[[paste0("SFC_", m)]] <- split_fusion_result2
}
    
    
one_result_dis <- do.call(rbind, all_results_dis)
rownames(one_result_dis) <- names(all_results_dis)
cat("Table: Ψ values\n")
print(round(one_result_dis, 3))
}

