function plot_NMI(NMI_MVKKM, NMI_MVSpec,NMI_MVKKM_paper, NMI_MVSpec_paper, p_list)

% Authors: Aaron Young, Alan Huang, Federico Berlfein
% April 23 2020
% COSI 126a Final Project

%function plots NMI results for comparison between our implementation and
%the paper's.
figure
hold on;

xv = 1:length(p_list);
plot(xv, NMI_MVKKM, 'r-x', xv, NMI_MVSpec,'b-*', xv, NMI_MVKKM_paper, 'k-+', xv, NMI_MVSpec_paper,'m-o'), legend('Our MVKKM', 'Our MVSpec', 'Paper MVKKM', 'Paper MVSpec'), xlabel('p value'),ylabel('NMI')%,axis([1 2 0 1.1])
xticks(xv)
legend('Location','southeast')
% first set of tick labels is for Digits data set. Second is for synthetic data
%xticklabels({'1','1.5','2','2.5','3','3.5', '4', '4.5', '5', '5.5', '6','Uniform'})
xticklabels({'1','1.3','1.5','2','4','6','Uniform'})
ylim([0.82 0.96])




hold off;