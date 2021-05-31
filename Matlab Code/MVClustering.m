function [Cluster_elem,w,Clustering_error]=MVClustering(Clusters, View_num, p, K, type)
%This code was based entirely from the Demo provided by G. Tzortzis
% We have simply made their demo into a function where we can run the
% "Demo" algorithm for any kernel, any p value, and any number of clusters.
% However, the actual code inside this function was based entirely on the
% one from G. Tzortzis, so all credit goes to him. We simply made a few
% adjustments to make it into a function.
%Courtesy of G. Tzortzis



%Algorithm parameters.
%---------------------
Init_w=ones(1,View_num)/View_num; %initial view weights (uniform as in the paper). 
%---------------------


%Cluster the instances using the MVKKM procedure.
%------------------------------------------------norma

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
if strcmp(type,'MVKKM')
    fprintf('Clustering with MVKKM for p=%f\n',p);
    fprintf('------------------------------------\n\n');
    [Cluster_elem,w,Clustering_error]=MVKKM(K,Clusters,View_num,p,Init_cluster_elem,Init_w);

    fprintf('Final MVKKM objective is:%f\n\n',Clustering_error);

    fprintf('Final view weights\n');
    for i=1:View_num
        fprintf('View %d: %f   ',i,w(i));
    end
    fprintf('\n\n');

    fprintf('End of clustering with MVKKM\n');
    fprintf('----------------------------\n\n');
    return
end
%------------------------------------------------

if strcmp(type,'MVSpec')
        %Cluster the instances using the MVSpec procedure.
    %------------------------------------------------

    fprintf('Clustering with MVSpec for p=%f\n',p);
    fprintf('-------------------------------------\n\n');

    %Perform multi-view clustering using MVSpec.
    %Initial clusters are not required for MVSpec.
    [Cluster_elem,w,Clustering_error]=MVSpec(K,Clusters,View_num,p,Init_w);

    %Plot the MVSpec solution.

    fprintf('Final MVSpec objective is:%f\n\n',Clustering_error);

    fprintf('Final view weights\n');
    for i=1:View_num
        fprintf('View %d: %f   ',i,w(i));
    end
    fprintf('\n\n');

    fprintf('End of clustering with MVSpec\n');
    fprintf('-----------------------------\n\n');

    %------------------------------------------------
    return
end
return