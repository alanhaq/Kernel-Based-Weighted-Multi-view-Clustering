function [Cluster_elem,w,Clustering_error]=MVSpec(K,Clusters,View_num,p,Init_w)
%
%[Cluster_elem,w,Clustering_error]=MVSpec(K,Clusters,View_num,p)
%
%This function implements the multi-view spectral clustering (MVSpec) algorithm as described in  
%G.Tzortzis and A.Likas, "Kernel-based Weighted Multi-view Clustering", ICDM 2012.
%
%Note that before applying MVSpec it is important to normalize the views to
%have comparable intra-cluster variances (refer to the paper for more details).
%
%Function Inputs
%===============
%
%K is a matrix containing the views' (positive-definite) kernel matrices.
%If V views are available and N data points then K is of size VNxN,
%i.e. the kernel matrices are stacked one below the other to form K.
%(K=[Kernel 1; Kernel 2;...;Kernel V]).
%
%Clusters is the number of clusters we search for.
%
%View_num is the number of views present in the dataset.
%
%p is the exponent that controls the sparsity of the kernel (view) weights.
%Smaller values lead to sparser solutions. Note that p>=1.
%
%Inti_w is a row vector containing initial values for the kernel (view)
%weights w_v (w_v>=0 & sum_{v=1}^V w_v = 1). A uniform initialization is
%used in the paper, but this can be changed if prior knowledge is available for the views.
%
%Function Outputs
%================
%
%Cluster_elem is a column vector containing the final discrete partitioning of the dataset. The clusters are indexed 1,...,Clusters.
%
%w is a row vector containing the final kernel (view) weights.
%
%Clustering_error is the value of the MVSpec objective function (relaxed intra-cluster variance).
%
%
%Courtesy of G. Tzortzis

if p<1
    error('p must be greater or equal to 1');
end

if sum(Init_w<0)~=0 || sum(Init_w)~=1
       error('Weights must be positive and sum to unity');
end

%Number of points in each view.
View_data_num=size(K,1)/View_num;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the implementation of MVSpec, we use the relaxed trace-based formulation of the intra-cluster
%variance objective as described by (8) in our paper.
%Our target is to minimize  sum_{v=1}^V w_v^p (tr(K_v) - tr(Y' K_v Y)), wrt Y and w_v  (sum_{v=1}^V w_v = 1).
%We alternate between optimizing Y given w_v and w_v given Y.
%For MVSpec Y is an arbitrary orthonormal matrix (relaxed cluster indicator matrix).


%Store the trace values tr(K_v) for each kernel K_v which are fixed and do not change.
Trace_per_kernel_1=zeros(View_num,1);
for i=1:View_num
    Trace_per_kernel_1(i)=trace(K((i-1)*View_data_num+1:i*View_data_num,:));
end

%Start iterations.
Iter=1;
Old_Clustering_error=inf;
w=Init_w;

while 1
    
    fprintf('+++++++++MVSpec Iteration %d+++++++++\n\n',Iter);
   
    %1st part - Update the clusters for given weights.
    %-------------------------------------------------
    fprintf('Updating the clusters\n');

    %Calculate the composite kernel (see (4) in the paper).
    K_sum=zeros(View_data_num);
    for i=1:View_num
        K_sum=K_sum+w(i)^p*K((i-1)*View_data_num+1:i*View_data_num,:);
    end

    %Calculate the top-"Clusters" eigenvectors of the composite kernel to get the relaxed cluster indicator matrix.
    %We do not search for a discrete partitioning here.
    [Y,eig_val]=eigs(K_sum,Clusters,'LA');
            
    %Store the new trace values tr(Y' K_v Y) for each kernel K_v after getting the relaxed partitioning.
    Trace_per_kernel=zeros(View_num,1);
    for i=1:View_num
        Trace_per_kernel(i)=trace(Y'*K((i-1)*View_data_num+1:i*View_data_num,:)*Y);
    end
    
    %Store the per-view trace differences.
    %These are the per-view relaxed intra-cluster variances after updating the clusters (D_v quantity in the paper).
    Trace_per_kernel_diff=Trace_per_kernel_1-Trace_per_kernel;
    
    %Calculate the objective (relaxed intra-cluster variance) value after updating the clusters.
    Clustering_error=(w.^p)*Trace_per_kernel_diff;
    fprintf('The objective after updating the clusters is:%f\n',Clustering_error);
    
    %Check for convergence.
    if abs(1-(Clustering_error/Old_Clustering_error))<0.0001
        fprintf('MVSpec reached convergence\n');
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
    fprintf('\n++++++++++++++++++++++++++++++++++++\n\n');
        
    Iter=Iter+1;
end

fprintf('\n++++++++++++++++++++++++++++++++++++\n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the final discrete partitioning.
%Cluster_elem contains the discrete cluster assignments.
Cluster_elem=Spectral_clustering(Clusters,K_sum);
 
return
