function plot_data(Clusters,Dataset,Cluster_elem)
%
%plot_data(Clusters,Dataset,Cluster_elem)
%
%This function plots 2-d points with different colors denoting the
%different clusters.
%
%Function Inputs
%===============
%
%Clusters is the number of clusters.
%
%Dataset is a matrix containing the 2-d instances. Each row corresponds to
%an instance.
%
%Cluster_elem is a vector containing the cluster indices. The clusters are 
%indexed 1,2,...,Clusters.
%
%
%Courtesy of G. Tzortzis

figure
hold on;

for i=1:Clusters
    This_elem=find(Cluster_elem==i);
    plot(Dataset(This_elem,1),Dataset(This_elem,2),'+','Color',[i/Clusters,1-i/Clusters,0]);     
end
axis equal;

hold off;

return