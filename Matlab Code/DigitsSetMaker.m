% Authors: Aaron Young, Alan Huang, Federico Berlfein
% April 23 2020
% COSI 126a Final Project

% Code produces the kernel's for all views in the Multiple Feature
% Handwriten Digits Data set presented in the paper. The 5 views are
% loaded, and 5 subsets of the data are created to emulate the paper's data
% sets. Each atribute column for every view is normalized to unit variance
% and views are then passed on through a gaussian image-filter that helps
% reduce noise. RBF kernel matrices are made for every view, and are then
% normalized based on the average pairwise square distance for the kernel
% For every subset of the data, a Kernel with each view's kernel is saved,
% where view's kernels are stacked on top of each other.

% download 5 different features from tables
%clear all

gamma = 0.001;
Clusters=4; %number of clusters.
View_num=5; %number of views present in the dataset.

View_prof_corr = table2array(readtable('mfeat-fac')); % profile correlations
View_fourier = table2array(readtable('mfeat-fou')); % fourier coeff
View_KL = table2array(readtable('mfeat-kar')); %Karhunen-Love coefficients
View_pix_ave = table2array(readtable('mfeat-kar')); % pixel averages
View_Z_mom = table2array(readtable('mfeat-kar')); % Zernike moments

% The first three lines of every view are the subset selection, the
% normalization to unit variance for every column, and the gaussian image
% filtering
% the last two lines for every view are the RBF kernel matrix calculation
% and normalization

% MF0169
prof_corr_0169 = [View_prof_corr(1:400,:); View_prof_corr(1201:1400, :);View_prof_corr(1801:2000, :) ];
prof_corr_0169 =  (prof_corr_0169) ./ repmat(std(prof_corr_0169), size(prof_corr_0169,1), 1);
prof_corr_0169 = imgaussfilt(prof_corr_0169); % gaussian image-filter
K_prof_corr_0169 = pairwise_kernels(prof_corr_0169, prof_corr_0169,  'rbf', gamma);
K_prof_corr_0169 = K_prof_corr_0169./mean(pdist(K_prof_corr_0169).^2);

fourier_0169 = [View_fourier(1:400,:); View_fourier(1201:1400, :);View_fourier(1801:2000, :) ];
fourier_0169 =  (fourier_0169) ./ repmat(std(fourier_0169), size(fourier_0169,1), 1);
fourier_0169 = imgaussfilt(fourier_0169);
K_fourier_0169 = pairwise_kernels(fourier_0169, fourier_0169,  'rbf', gamma);
K_fourier_0169 = K_fourier_0169./mean(pdist(K_fourier_0169).^2);

KL_0169 = [View_KL(1:400,:); View_KL(1201:1400, :);View_KL(1801:2000, :) ];
KL_0169 =  (KL_0169) ./ repmat(std(KL_0169), size(KL_0169,1), 1);
KL_0169 = imgaussfilt(KL_0169);
K_KL_0169 = pairwise_kernels(KL_0169, KL_0169,  'rbf', gamma);
K_KL_0169 = K_KL_0169./mean(pdist(K_KL_0169).^2);

pix_ave_0169 = [View_pix_ave(1:400,:); View_pix_ave(1201:1400, :);View_pix_ave(1801:2000, :) ];
pix_ave_0169 =  (pix_ave_0169) ./ repmat(std(pix_ave_0169), size(pix_ave_0169,1), 1);
pix_ave_0169 = imgaussfilt(pix_ave_0169);
K_pix_ave_0169 = pairwise_kernels(pix_ave_0169, pix_ave_0169,  'rbf', gamma);
K_pix_ave_0169 = K_pix_ave_0169./mean(pdist(K_pix_ave_0169).^2);

Z_mom_0169 = [View_Z_mom(1:400,:); View_Z_mom(1201:1400, :);View_Z_mom(1801:2000, :) ];
Z_mom_0169 =  (Z_mom_0169) ./ repmat(std(Z_mom_0169), size(Z_mom_0169,1), 1);
Z_mom_0169 = imgaussfilt(Z_mom_0169);
K_Z_mom_0169 = pairwise_kernels(Z_mom_0169, Z_mom_0169,  'rbf', gamma);
K_Z_mom_0169 = K_Z_mom_0169./mean(pdist(K_Z_mom_0169).^2);

% We now save the views and the Kernels separetely.
MF0169_rbf_kernel = [K_prof_corr_0169; K_fourier_0169;K_KL_0169;K_pix_ave_0169;K_Z_mom_0169 ];
save('MF0169.mat','prof_corr_0169', 'fourier_0169','K_KL_0169' ,'pix_ave_0169' ,'Z_mom_0169');
save('MF0169_rbf_kernel.mat', 'MF0169_rbf_kernel');

% this is repeated in the same manner for all subsets.

% MF0236
prof_corr_0236 = [View_prof_corr(1:200,:); View_prof_corr(401:800, :);View_prof_corr(1201:1400, :) ];
prof_corr_0236 =  (prof_corr_0236) ./ repmat(std(prof_corr_0236), size(prof_corr_0236,1), 1);
prof_corr_0236 = imgaussfilt(prof_corr_0236);
K_prof_corr_0236 = pairwise_kernels(prof_corr_0236, prof_corr_0236,  'rbf', gamma);
K_prof_corr_0236 = K_prof_corr_0236./mean(pdist(K_prof_corr_0236).^2);

fourier_0236 = [View_fourier(1:200,:); View_fourier(401:800, :);View_fourier(1201:1400, :) ];
fourier_0236 =  (fourier_0236) ./ repmat(std(fourier_0236), size(fourier_0236,1), 1);
fourier_0236 = imgaussfilt(fourier_0236);
K_fourier_0236 = pairwise_kernels(fourier_0236, fourier_0236,  'rbf', gamma);
K_fourier_0169 = K_fourier_0169./mean(pdist(K_fourier_0169).^2);

KL_0236 = [View_KL(1:200,:); View_KL(401:800, :);View_KL(1201:1400, :) ];
KL_0236 =  (KL_0236) ./ repmat(std(KL_0236), size(KL_0236,1), 1);
KL_0236 = imgaussfilt(KL_0236);
K_KL_0236 = pairwise_kernels(KL_0236, KL_0236,  'rbf', gamma);
K_KL_0236 = K_KL_0236./mean(pdist(K_KL_0236).^2);

pix_ave_0236 = [View_pix_ave(1:200,:); View_pix_ave(401:800, :);View_pix_ave(1201:1400, :) ];
pix_ave_0236 =  (pix_ave_0236) ./ repmat(std(pix_ave_0236), size(pix_ave_0236,1), 1);
pix_ave_0236 = imgaussfilt(pix_ave_0236);
K_pix_ave_0236 = pairwise_kernels(pix_ave_0236, pix_ave_0236,  'rbf', gamma);
K_pix_ave_0236 = K_pix_ave_0236./mean(pdist(K_pix_ave_0236).^2);

Z_mom_0236 = [View_Z_mom(1:200,:); View_Z_mom(401:800, :);View_Z_mom(1201:1400, :) ];
Z_mom_0236 =  (Z_mom_0236) ./ repmat(std(Z_mom_0236), size(Z_mom_0236,1), 1);
Z_mom_0236 = imgaussfilt(Z_mom_0236);
K_Z_mom_0236 = pairwise_kernels(Z_mom_0236, Z_mom_0236,  'rbf', gamma);
K_Z_mom_0236 = K_Z_mom_0236./mean(pdist(K_Z_mom_0236).^2);


MF0236_rbf_kernel = [K_prof_corr_0236; K_fourier_0236;K_KL_0236;K_pix_ave_0236;K_Z_mom_0236 ];
save('MF0236.mat','prof_corr_0236', 'fourier_0236','K_KL_0236' ,'pix_ave_0236' ,'Z_mom_0236');
save('MF0236_rbf_kernel.mat', 'MF0236_rbf_kernel');

% MF1456
prof_corr_1456 = [View_prof_corr(201:400,:); View_prof_corr(801:1400, :) ];
prof_corr_1456 =  (prof_corr_1456) ./ repmat(std(prof_corr_1456), size(prof_corr_1456,1), 1);
prof_corr_1456 = imgaussfilt(prof_corr_1456);
K_prof_corr_1456 = pairwise_kernels(prof_corr_1456, prof_corr_1456,  'rbf', gamma);
K_prof_corr_1456 = K_prof_corr_1456./mean(pdist(K_prof_corr_1456).^2);

fourier_1456 = [View_fourier(201:400,:); View_fourier(801:1400, :) ];
fourier_1456 =  (fourier_1456) ./ repmat(std(fourier_1456), size(fourier_1456,1), 1);
fourier_1456 = imgaussfilt(fourier_1456);
K_fourier_1456 = pairwise_kernels(fourier_1456, fourier_1456,  'rbf', gamma);
K_fourier_1456 = K_fourier_1456./mean(pdist(K_fourier_1456).^2);

KL_1456 = [View_KL(201:400,:); View_KL(801:1400, :) ];
KL_1456 =  (KL_1456) ./ repmat(std(KL_1456), size(KL_1456,1), 1);
KL_1456 = imgaussfilt(KL_1456);
K_KL_1456 = pairwise_kernels(KL_1456, KL_1456,  'rbf', gamma);
K_KL_1456 = K_KL_1456./mean(pdist(K_KL_1456).^2);

pix_ave_1456 = [View_pix_ave(201:400,:); View_pix_ave(801:1400, :) ];
pix_ave_1456 =  (pix_ave_1456) ./ repmat(std(pix_ave_1456), size(pix_ave_1456,1), 1);
pix_ave_1456 = imgaussfilt(pix_ave_1456);
K_pix_ave_1456 = pairwise_kernels(pix_ave_1456, pix_ave_1456,  'rbf', gamma);
K_pix_ave_1456 = K_pix_ave_1456./mean(pdist(K_pix_ave_1456).^2);

Z_mom_1456 = [View_Z_mom(201:400,:); View_Z_mom(801:1400, :)];
Z_mom_1456 =  (Z_mom_1456) ./ repmat(std(Z_mom_1456), size(Z_mom_1456,1), 1);
Z_mom_1456 = imgaussfilt(Z_mom_1456);
K_Z_mom_1456 = pairwise_kernels(Z_mom_1456, Z_mom_1456,  'rbf', gamma);
K_Z_mom_1456 = K_Z_mom_1456./mean(pdist(K_Z_mom_1456).^2);


MF1456_rbf_kernel = [K_prof_corr_1456; K_fourier_1456;K_KL_1456;K_pix_ave_1456;K_Z_mom_1456 ];
save('MF1456.mat','prof_corr_1456', 'fourier_1456','K_KL_1456' ,'pix_ave_1456' ,'Z_mom_1456');
save('MF1456_rbf_kernel.mat', 'MF1456_rbf_kernel');


% MF1367
prof_corr_1367 = [View_prof_corr(201:400,:); View_prof_corr(601:800, :);View_prof_corr(1201:1600, :) ];
prof_corr_1367 =  (prof_corr_1367) ./ repmat(std(prof_corr_1367), size(prof_corr_1367,1), 1);
prof_corr_1367 = imgaussfilt(prof_corr_1367);
K_prof_corr_1367 = pairwise_kernels(prof_corr_1367, prof_corr_1367,  'rbf', gamma);
K_prof_corr_1367 = K_prof_corr_1367./mean(pdist(K_prof_corr_1367).^2);

fourier_1367 = [View_fourier(201:400,:); View_fourier(601:800, :);View_fourier(1201:1600, :) ];
fourier_1367 =  (fourier_1367) ./ repmat(std(fourier_1367), size(fourier_1367,1), 1);
fourier_1367 = imgaussfilt(fourier_1367);
K_fourier_1367 = pairwise_kernels(fourier_1367, fourier_1367,  'rbf', gamma);
K_fourier_1367 = K_fourier_1367./mean(pdist(K_fourier_1367).^2);

KL_1367 = [View_KL(201:400,:); View_KL(601:800, :);View_KL(1201:1600, :) ];
KL_1367 =  (KL_1367) ./ repmat(std(KL_1367), size(KL_1367,1), 1);
KL_1367 = imgaussfilt(KL_1367);
K_KL_1367 = pairwise_kernels(KL_1367, KL_1367,  'rbf', gamma);
K_KL_1367 = K_KL_1367./mean(pdist(K_KL_1367).^2);

pix_ave_1367 = [View_pix_ave(201:400,:); View_pix_ave(601:800, :);View_pix_ave(1201:1600, :) ];
pix_ave_1367 =  (pix_ave_1367) ./ repmat(std(pix_ave_1367), size(pix_ave_1367,1), 1);
pix_ave_1367 = imgaussfilt(pix_ave_1367);
K_pix_ave_1367 = pairwise_kernels(pix_ave_1367, pix_ave_1367,  'rbf', gamma);
K_pix_ave_1367 = K_pix_ave_1367./mean(pdist(K_pix_ave_1367).^2);

Z_mom_1367 = [View_Z_mom(201:400,:); View_Z_mom(601:800, :);View_Z_mom(1201:1600, :) ];
Z_mom_1367 =  (Z_mom_1367) ./ repmat(std(Z_mom_1367), size(Z_mom_1367,1), 1);
Z_mom_1367 = imgaussfilt(Z_mom_1367);
K_Z_mom_1367 = pairwise_kernels(Z_mom_1367, Z_mom_1367,  'rbf', gamma);
K_Z_mom_1367 = K_Z_mom_1367./mean(pdist(K_Z_mom_1367).^2);


MF1367_rbf_kernel = [K_prof_corr_1367; K_fourier_1367; K_KL_1367; K_pix_ave_1367;K_Z_mom_1367 ];
save('MF1367.mat','prof_corr_1367', 'fourier_1367','K_KL_1367' ,'pix_ave_1367' ,'Z_mom_1367');
save('MF1367_rbf_kernel.mat', 'MF1367_rbf_kernel');


% MF4689
prof_corr_4689 = [View_prof_corr(801:1000,:); View_prof_corr(1201:1400, :);View_prof_corr(1601:2000, :) ];
prof_corr_4689 =  (prof_corr_4689) ./ repmat(std(prof_corr_4689), size(prof_corr_4689,1), 1);
prof_corr_4689 = imgaussfilt(prof_corr_4689);
K_prof_corr_4689 = pairwise_kernels(prof_corr_4689, prof_corr_4689,  'rbf', gamma);
K_prof_corr_4689 = K_prof_corr_4689./mean(pdist(K_prof_corr_4689).^2);

fourier_4689 = [View_fourier(801:1000,:); View_fourier(1201:1400, :);View_fourier(1601:2000, :) ];
fourier_4689 =  (fourier_4689) ./ repmat(std(fourier_4689), size(fourier_4689,1), 1);
fourier_4689 = imgaussfilt(fourier_4689);
K_fourier_4689 = pairwise_kernels(fourier_4689, fourier_4689,  'rbf', gamma);
K_fourier_4689 = K_fourier_4689./mean(pdist(K_fourier_4689).^2);

KL_4689 = [View_KL(801:1000,:); View_KL(1201:1400, :);View_KL(1601:2000, :) ];
KL_4689 =  (KL_4689) ./ repmat(std(KL_4689), size(KL_4689,1), 1);
KL_4689 = imgaussfilt(KL_4689);
K_KL_4689 = pairwise_kernels(KL_4689, KL_4689,  'rbf', gamma);
K_KL_4689 = K_KL_4689./mean(pdist(K_KL_4689).^2);

pix_ave_4689 = [View_pix_ave(801:1000,:); View_pix_ave(1201:1400, :);View_pix_ave(1601:2000, :) ];
pix_ave_4689 =  (pix_ave_4689) ./ repmat(std(pix_ave_4689), size(pix_ave_4689,1), 1);
pix_ave_4689 = imgaussfilt(pix_ave_4689);
K_pix_ave_4689 = pairwise_kernels(pix_ave_4689, pix_ave_4689,  'rbf', gamma);
K_pix_ave_4689 = K_pix_ave_4689./mean(pdist(K_pix_ave_4689).^2);

Z_mom_4689 = [View_Z_mom(801:1000,:); View_Z_mom(1201:1400, :);View_Z_mom(1601:2000, :) ];
Z_mom_4689 =  (Z_mom_4689) ./ repmat(std(Z_mom_4689), size(Z_mom_4689,1), 1);
Z_mom_4689 = imgaussfilt(Z_mom_4689);
K_Z_mom_4689 = pairwise_kernels(Z_mom_4689, Z_mom_4689,  'rbf', gamma);
K_Z_mom_4689 = K_Z_mom_4689./mean(pdist(K_Z_mom_4689).^2);


MF4689_rbf_kernel = [K_prof_corr_4689; K_fourier_4689; K_KL_4689; K_pix_ave_4689;K_Z_mom_4689 ];
save('MF4689.mat','prof_corr_4689', 'fourier_4689','K_KL_4689' ,'pix_ave_4689' ,'Z_mom_4689');
save('MF4689_rbf_kernel.mat', 'MF4689_rbf_kernel');

