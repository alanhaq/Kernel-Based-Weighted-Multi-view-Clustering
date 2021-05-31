% Authors: Aaron Young, Alan Huang, Federico Berlfein
% April 23 2020
% COSI 126a Final Project

%Code loads and plots MVKKM and MVSpec results for Paper's Synthetic
% dataset. This data set is already pre-computed by the paper's authors and
% provides an easy check on the functionality of the code.


load synthetic_data.mat;
%load synthetic_data_kernel.mat;
gamma = 12.5;

%View_1 =  (View_1) ./ repmat(std(View_1), size(View_1,1), 1);
%View_1 = imgaussfilt(View_1); % gaussian image-filter
K_View_1 = pairwise_kernels(View_1, View_1,  'rbf', gamma);
%K_View_1 = K_View_1./mean(pdist(K_View_1).^2); 

%View_2 =  (View_2) ./ repmat(std(View_2), size(View_2,1), 1);
%View_2 = imgaussfilt(View_2); % gaussian image-filter
K_View_2 = pairwise_kernels(View_2, View_2,  'rbf', gamma);
%K_View_2 = K_View_2./mean(pdist(K_View_2).^2);

K = [K_View_1;K_View_2];

p_list = [1 1.3 1.5 2 4 6 30];
gamma = 0.2;


Clusters=3; %number of clusters.
View_num=2; %number of views present in the dataset.

Cluster_labels_MVKKM = [];
w_list_MVKKM = [];
Clustering_errors_MVKKM = [];
NMI_MVKKM = [];

Cluster_labels_MVSpec = [];
w_list_MVSpec = [];
Clustering_errors_MVSpec = [];
NMI_MVSpec = [];

% test different values of p for both algorithms
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

% display NMI results for both algorithms and final weights
disp('NMI_MVKKM:  ');
disp(NMI_MVKKM);
disp('NMI_MVSpec:  ');
disp(NMI_MVSpec);
disp('MVKKM Weights:  ');
disp(w_list_MVKKM);
disp('MVSpec Weights:  ');
disp(w_list_MVSpec);

NMI_MVKKM_paper = [1 1 1 0.769 0.749 0.747 0.701 ];
NMI_MVSpec_paper = [0.681 0.671 0.663 0.632 0.593 0.593 0.552];
plot_NMI(NMI_MVKKM,NMI_MVSpec, NMI_MVKKM_paper, NMI_MVSpec_paper, p_list) % plot our results vs papers
