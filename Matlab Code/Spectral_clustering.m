function [Cluster_elem]=Spectral_clustering(Clusters,K)
%
%[Cluster_elem]=Spectral_clustering(Clusters,K)
%
%This function performs spectral clustering on the affinity matrix K.
%In the case of the MVSpec algorithm the affinity matrix is the composite kernel of the views.
%The implementation follows the steps of the spectral clustering algorithm described in
%the paper: A. Ng, M. Jordan and Y. Weiss, "On Spectral Clustering: Analysis and an Algorithm", NIPS, 2001.
%The only difference is that the eigen-analysis is directly performed on
%the affinity matrix, instead of its Laplacian as in Ng et al.
%This is done as according to the MVSpec method the eigenvectors of the
%composite kernel are needed.
%
%Function Inputs
%===============
%
%Clusters is the number of clusters we search for.
%
%K is the dataset (square) affinity matrix. For the MVSpec algorithm the
%affinity matrix is the composite kernel of the views.
%
%Function Outputs
%================
%
%Cluster_elem is a column vector containing the final discrete partitioning of the dataset. The clusters are indexed 1,...,Clusters.
%
%
%Courtesy of G. Tzortzis

rand('state',0);

%Calculate the top-"Clusters" eigenvectors to form the relaxed cluster indicator matrix.
%These eigenvectors are an embedding of the original instances in a "Clusters"-dimensional space.
[eig_vec,eig_val]=eigs(K,Clusters,'LA');

%Normalize the embedding of each instance to unit length.
normalize=sqrt(sum(eig_vec.^2,2));
eig_vec=eig_vec./repmat(normalize,1,Clusters);

%Cluster the embeddings using k-means to get a discrete partitioning of the data.
%K-means is restarted 30 times.
Cluster_elem=kmeans(eig_vec,Clusters,'Distance','sqEuclidean','Start','sample','Maxiter',1000,'EmptyAction','singleton','Display','off','Replicates',30);

return



