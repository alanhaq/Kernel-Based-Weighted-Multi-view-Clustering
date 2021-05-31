% Authors: Aaron Young, Alan Huang, Federico Berlfein
% April 23 2020
% COSI 126a Final Project

% Code loads Kernels for views of Multiple Feautre Written Digits Data set
% for subset MF0169. Similarly to SyntheticData.m, this caluclates the NMI
% for different p values for MVKKM and MVSpec, and compares our results to
% the ones from the paper.


load MF0169_rbf_kernel.mat;

K = MF0169_poly_kernel;
%gamma = 0.2;
Clusters=4; %number of clusters.
View_num=5; %number of views present in the dataset.
View_data_num=size(K,1)/View_num;

Ground_truth = [ones(200,1); 2.* ones(200,1); 3.* ones(200,1); 4.* ones(200,1) ];

p_list = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 50];

Cluster_labels_MVKKM = [];
w_list_MVKKM = [];
Clustering_errors_MVKKM = [];
NMI_MVKKM = [];

Cluster_labels_MVSpec = [];
w_list_MVSpec = [];
Clustering_errors_MVSpec = [];
NMI_MVSpec = [];

for i=1:length(p_list(:))
    %MVKKM
    [Cluster_elem_MVKKM,w_MVKKM,Clustering_error_MVKKM] = MVClustering(Clusters, View_num, p_list(i), K, 'MVKKM');
    Cluster_labels_MVKKM = [Cluster_labels_MVKKM Cluster_elem_MVKKM]; % colms represent dif values of p
    w_list_MVKKM = [w_list_MVKKM ;w_MVKKM]; % rows represent dif values of p
    Clustering_errors_MVKKM = [Clustering_errors_MVKKM Clustering_error_MVKKM];  % colms represent dif values of p
    nmi_elem_MVKKM = nmi(Cluster_elem_MVKKM, Ground_truth);
    NMI_MVKKM = [NMI_MVKKM nmi_elem_MVKKM];
    
    %MVSpec
    [Cluster_elem_MVSpec,w_MVSpec,Clustering_error_MVSpec] = MVClustering(Clusters, View_num, p_list(i), K, 'MVSpec');
    Cluster_labels_MVSpec = [Cluster_labels_MVSpec Cluster_elem_MVSpec]; % colms represent dif values of p
    w_list_MVSpec = [w_list_MVSpec ;w_MVSpec]; % rows represent dif values of p
    Clustering_errors_MVSpec = [Clustering_errors_MVSpec Clustering_error_MVSpec];  % colms represent dif values of p
    nmi_elem_MVSpec = nmi(Cluster_elem_MVSpec, Ground_truth);
    NMI_MVSpec = [NMI_MVSpec nmi_elem_MVSpec];
end

% display results
disp('NMI_MVKKM:  ');
disp(NMI_MVKKM);
disp('NMI_MVSpec:  ');
disp(NMI_MVSpec);
disp('MVKKM Weights:  ');
disp(w_list_MVKKM);
disp('MVSpec Weights:  ');
disp(w_list_MVSpec);

%plot results vs paper's
NMI_MVKKM_paper = [0.905 0.94 0.95 0.94 0.93 0.927 0.925 0.925 0.925 0.925 0.925 0.92];
NMI_MVSpec_paper = [0.84 0.88 0.90 0.895 0.895 0.895 0.89 0.89 0.89 0.89 0.89 0.885];
plot_NMI(NMI_MVKKM,NMI_MVSpec,NMI_MVKKM_paper,NMI_MVSpec_paper, p_list)
%title('MVKKM Solution');
