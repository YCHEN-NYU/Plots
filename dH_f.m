
clear;
close all;

% Define color table & marker table below for plotting
colortable = ['r','b','k','c','g','m','r','b',...
    'k','c','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^',...
    'o','s','v','^','o','s','v','^'];
sampleName = 'STT5444';

titlename = ['$f - H(' sampleName  ')$'];
xlabelname = '$H_{res} (T)$';
ylabelname = '$f (GHz)$';


yleft = 7;
yright = 27.5;


HLowerbound = 0;
HUpperbound = 1;

fLowerbound = 0;
fUpperbound = 30;

dHLowerbound = 0;
dHUpperbound = 150;

Hlim = [HLowerbound,HUpperbound]; % Hres(T)
flim = [fLowerbound,fUpperbound]; % f(GHz)
dHlim = [dHLowerbound, dHUpperbound];
f2lim = [0, 10];

fileFormat = ['*', sampleName '*.txt'];
files=dir(fileFormat);
[filenames, index] = sort_nat({files.name});

fig = figure();
set(fig, 'Position', [200, 100, 1000, 800]);
set(fig,'color','w');
fig.PaperPositionMode = 'auto';% set image size as auto
figfH = 'fH.png';
N = 1;

plot1 = zeros(1,2*N);
for i = 1:2*N
    rawdata = importdata(filenames{i});
    
    f = rawdata(:,1);
    H = rawdata(:,2)/10^4;
    lw = rawdata(:,3);
    lw_lb = rawdata(:,4);
    lw_ub = rawdata(:,5);
    
%     if i > 3 
%         H = H+0.032; % calibration of Hall sensor difference
%     end

    x = H;
    y = f;
%     y_err = (lw_ub -lw_lb)/2;
    if i <=N
    plot1(i) = plot(x,y,'color',colortable(i), ...
    'marker',markertable(i),'markersize',10);
    else 
        plot1(i) = plot(x,y,'color',colortable(i), ...
    'marker',markertable(i),'markersize',10,'MarkerFaceColor',colortable(i));
    end
hold on;

title(titlename,'FontSize',42,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlim(Hlim);
ylim(flim);


ylabel(ylabelname,'FontSize',36,'FontWeight',...
    'bold', 'interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlabel(xlabelname,'FontSize',36,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',36,'Linewidth',3,'fontweight','bold');

end


figure(fig);
% legend(plot1, '2.67 nm - old', '2.06 nm - old','1.69 nm - old','2.67 nm - new', '2.06 nm - new','1.69 nm - new', ...
%      'location', 'east');

legend(plot1, '4.0 nm - old', '4.0 nm - new','location', 'east');


 fignameKittel = [sampleName '_Kittel.png'];

 saveas(fig,fignameKittel);
