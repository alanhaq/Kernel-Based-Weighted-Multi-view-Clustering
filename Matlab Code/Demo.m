%This demo shows how to call the MVKKM and MVSpec functions that implement
%the multi-view kernel k-means and multi-view spectral clustering methods,
%respectively, described in the paper: G.Tzortzis and A.Likas, 
%"Kernel-based Weighted Multi-view Clustering", ICDM 2012.
%For the demonstration, the synthetic dataset of the above paper is used.
%It consists of two views (the second view is a noisy version of the first)
%and three clusters. A precomputed rbf kernel is utilized for each view
%(sigma=0.2 for both views), as in the experiments.
%
%
%Courtesy of G. Tzortzis

clear all
close all

%Load the instances.
%Matrices View_1 and View_2 contain the instances, matrix Ground_truth
%contains the ground truth for the dataset.
load synthetic_data.mat;

%Plot the two views.
plot_data(3,View_1,Ground_truth);
title('View 1');
plot_data(3,View_2,Ground_truth);
title('View 2 - Noisy View');
drawnow;


%Algorithm parameters.
%---------------------
Clusters=3; %number of clusters.
View_num=2; %number of views present in the dataset.
p=1.5; %the exponent p (p>=1) (smaller p values lead to sparser view weights).
Init_w=ones(1,View_num)/View_num; %initial view weights (uniform as in the paper). 
%---------------------


%Load the precomputed rbf kernel for each view.
%Both view kernels are contained in a single matrix, K, one below the other.
load synthetic_data_kernel.mat;


%Cluster the instances using the MVKKM procedure.
%------------------------------------------------norma

fprintf('Clustering with MVKKM for p=%f\n',p);
fprintf('------------------------------------\n\n');

View_data_num=size(K,1)/View_num;

%Calculate the composite kernel (see (4) in the paper).
%The composite kernel is used to get an initial partitioning with the
%(single-view) global kernel k-means algorithm.
K_sum=zeros(View_data_num);
for i=1:View_num
    K_sum=K_sum+Init_w(i)^p*K((i-1)*View_data_num+1:i*View_data_num,:);
end
   
fprintf('Global kernel k-means initialization\n');

%Before applying MVKKM, initial clusters must be found.
%Here the incremental-deterministic global kernel k-means algorithm is utilized for this purpose (as in the paper).
Init_cluster_elem=Weighted_Global_Kernel_K_Means(K_sum,ones(View_data_num,1),Clusters,'-');
fprintf('End of initialization\n\n');

%Perform multi-view clustering using MVKKM.
[Cluster_elem_MVKKM,w_MVKKM,Clustering_error_MVKKM]=MVKKM(K,Clusters,View_num,p,Init_cluster_elem,Init_w);

%Plot the MVKKM solution.
plot_data(3,View_1,Cluster_elem_MVKKM);
title('MVKKM Solution');

fprintf('Final MVKKM objective is:%f\n\n',Clustering_error_MVKKM);

fprintf('Final view weights\n');
for i=1:View_num
    fprintf('View %d: %f   ',i,w_MVKKM(i));
end
fprintf('\n\n');

fprintf('End of clustering with MVKKM\n');
fprintf('----------------------------\n\n');

%------------------------------------------------


%Cluster the instances using the MVSpec procedure.
%------------------------------------------------

fprintf('Clustering with MVSpec for p=%f\n',p);
fprintf('-------------------------------------\n\n');

%Perform multi-view clustering using MVSpec.
%Initial clusters are not required for MVSpec.
[Cluster_elem_MVSpec,w_MVSpec,Clustering_error_MVSpec]=MVSpec(K,Clusters,View_num,p,Init_w);

%Plot the MVSpec solution.
plot_data(3,View_1,Cluster_elem_MVSpec);
title('MVSpec Solution');

fprintf('Final MVSpec objective is:%f\n\n',Clustering_error_MVSpec);

fprintf('Final view weights\n');
for i=1:View_num
    fprintf('View %d: %f   ',i,w_MVSpec(i));
end
fprintf('\n\n');

fprintf('End of clustering with MVSpec\n');
fprintf('-----------------------------\n\n');

%------------------------------------------------
