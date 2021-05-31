function [Cluster_elem,w,Clustering_error]=MVKKM(K,Clusters,View_num,p,Init_cluster_elem,Init_w)
%
%[Cluster_elem,w,Clustering_error]=MVKKM(Clusters,View_num,p,Init_cluster_elem)
%
%This function implements the multi-view kernel k-means (MVKKM) algorithm as described in  
%G.Tzortzis and A.Likas, "Kernel-based Weighted Multi-view Clustering", ICDM 2012.
%
%Note that before applying MVKKM it is important to normalize the views to
%have comparable intra-cluster variances (refer to the paper for more details).
%
%Function Inputs
%===============
%
%K is a matrix containing the views' (positive-definite) kernel matrices.
%If V views are available and N data points then K is of size VNxN,
%i.e. the kernel matrices are stacked one below the other to form K 
%(K=[Kernel 1;Kernel 2;...;Kernel V]).
%
%Clusters is the number of clusters we search for.
%
%View_num is the number of views present in the dataset.
%
%p is the exponent that controls the sparsity of the kernel (view) weights.
%Smaller values lead to sparser solutions. Note that p>=1.
%
%Init_cluster_elem is an initial partitioning of the dataset. It is a column vector
%containing the cluster index of each point. The clusters are indexed
%1,2,...,Clusters. In the paper we applied global kernel k-means to get
%initial cluster assignments, but any other algorithm can be used instead.
%
%Inti_w is a row vector containing initial values for the kernel (view)
%weights w_v (w_v>=0 & sum_{v=1}^V w_v = 1). A uniform initialization is
%used in the paper, but this can be changed if prior knowledge is available for the views.
%
%Function Outputs
%================
%
%Cluster_elem is a column vector containing the final partitioning of the dataset. The clusters are indexed 1,...,Clusters.
%
%w is a row vector containing the final kernel (view) weights.
%
%Clustering_error is the value of the MVKKM objective function (intra-cluster variance).
%
%
%Courtesy of G. Tzortzis

if p<1
    error('p must be greater or equal to 1');
end

if sum(Init_w<0)~=0 || abs(sum(Init_w)-1)>10^(-15)
       error('Weights must be positive and sum to unity');
end

%Number of points in each view.
View_data_num=size(K,1)/View_num;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the implementation of MVKKM, we use the non-relaxed trace-based formulation 
%of the intra-cluster variance objective as described by (8) in our paper. 
%This formulation is equivalent to the distance-based formulation (5).
%Our target is to minimize  sum_{v=1}^V w_v^p (tr(K_v) - tr(Y' K_v Y)), wrt Y and w_v  (sum_{v=1}^V w_v = 1).
%We alternate between optimizing Y given w_v and w_v given Y.
%For MVKKM Y is an indicator matrix of the form described in (3).


%Store the trace values tr(K_v) for each kernel K_v which are fixed and do not change.
Trace_per_kernel_1=zeros(View_num,1);
for i=1:View_num
    Trace_per_kernel_1(i)=trace(K((i-1)*View_data_num+1:i*View_data_num,:));
end

%Start iterations.
Iter=1;
Old_Clustering_error=inf;
Cluster_elem=Init_cluster_elem;
w=Init_w;

while 1
    
    fprintf('+++++++++MVKKM Iteration %d+++++++++\n\n',Iter);
    
    %1st part - Update the clusters for given weights.
    %-------------------------------------------------
    fprintf('Updating the clusters\n');

    %Calculate the composite kernel (see (4) in the paper).
    K_sum=zeros(View_data_num);
    for i=1:View_num
        K_sum=K_sum+w(i)^p*K((i-1)*View_data_num+1:i*View_data_num,:);
    end

    %Run kernel k-means to get new assignments.
    %Kernel k-means is initialized using the partitioning returned by the previous iteration. 
    %Cluster_elem contains the new assignments.
    %Clustering_error is the objective (intra-cluster variance) value after updating the clusters.
    [Cluster_elem,Clustering_error]=Weighted_Kernel_K_Means(Cluster_elem,K_sum,ones(View_data_num,1),Clusters,'-');
    fprintf('The objective after updating the clusters is:%f\n',Clustering_error);
    
    if length(unique(Cluster_elem))<Clusters
        error('Empty clusters detected');
    end
    
    %Calculate the cluster indicator matrix Y (see (3) in the paper).
    Y=spalloc(View_data_num,Clusters,View_data_num);
    for i=1:Clusters
        Points_cluster_i=find(Cluster_elem==i);
        Y(Points_cluster_i,i)=1/sqrt(length(Points_cluster_i));
    end
    
    %Store the new trace values tr(Y' K_v Y) for each kernel K_v after the kernel k-means run.
    Trace_per_kernel=zeros(View_num,1);
    for i=1:View_num
        Trace_per_kernel(i)=trace(Y'*K((i-1)*View_data_num+1:i*View_data_num,:)*Y);
    end
    
    %Store the per-view trace differences.
    %These are the per-view intra-cluster variances after updating the clusters (D_v quantity in the paper).
    Trace_per_kernel_diff=Trace_per_kernel_1-Trace_per_kernel;
      
    %Check for convergence.
    if abs(1-(Clustering_error/Old_Clustering_error))<0.0001
        fprintf('MVKKM reached convergence\n');
        break
    end
    
    Old_Clustering_error=Clustering_error;

    %2nd part - Update the weights for given clusters.
    %------------------------------------------------
    fprintf('\nUpdating the weights\n');

    %Update the weights using the closed form expression (see (9) and (10) in the paper).
    if p~=1
        w=1./(Trace_per_kernel_diff.^(1/(p-1)).*sum(1./(Trace_per_kernel_diff.^(1/(p-1)))))';
    else
        [min_tr,min_tr_pos]=min(Trace_per_kernel_diff);
        w=zeros(1,View_num);
        w(min_tr_pos)=1;
    end
    
    %Tiny weight values are set to zero for better stability.
    w(w<10^-5)=0;
    w=w/sum(w);
            
    fprintf('The objective after updating the weights is:%f\n',(w.^p)*Trace_per_kernel_diff);
    fprintf('\n+++++++++++++++++++++++++++++++++++\n\n');
    
    Iter=Iter+1;
end

fprintf('\n+++++++++++++++++++++++++++++++++++\n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
