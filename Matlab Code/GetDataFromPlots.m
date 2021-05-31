%MF01690
open MF0169Overplot.fig
h = findobj(gca,'Type','line');
%x=get(h,'Xdata');
y=get(h,'Ydata');
NMI_MVKKM_0169 = y{3,1};
NMI_MVSpec_0169 = y{4,1};
save('MF0169_NMI.mat','NMI_MVKKM_0169', 'NMI_MVKKM_0169');

%MF0236
open MF0236Overplot.fig
h = findobj(gca,'Type','line');
%x=get(h,'Xdata');
y=get(h,'Ydata');
NMI_MVKKM_0236 = y{3,1};
NMI_MVSpec_0236 = y{4,1};
save('MF0236_NMI.mat','NMI_MVKKM_0236', 'NMI_MVKKM_0236');

%MF1367
open MF1367Overplot.fig
h = findobj(gca,'Type','line');
%x=get(h,'Xdata');
y=get(h,'Ydata');
NMI_MVKKM_1367 = y{3,1};
NMI_MVSpec_1367 = y{4,1};
save('MF1367_NMI.mat','NMI_MVKKM_1367', 'NMI_MVKKM_1367');

%MF1456
open MF1456Overplot.fig
h = findobj(gca,'Type','line');
%x=get(h,'Xdata');
y=get(h,'Ydata');
NMI_MVKKM_1456 = y{3,1};
NMI_MVSpec_1456 = y{4,1};
save('MF1456_NMI.mat','NMI_MVKKM_1456', 'NMI_MVKKM_1456');

%MF4689
open MF4689Overplot.fig
h = findobj(gca,'Type','line');
%x=get(h,'Xdata');
y=get(h,'Ydata');
NMI_MVKKM_4689 = y{3,1};
NMI_MVSpec_4689 = y{4,1};
save('MF4689_NMI.mat','NMI_MVKKM_4689', 'NMI_MVKKM_4689');
