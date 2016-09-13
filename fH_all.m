% ========================================
clear;
close all;
% ========================================
% Define color table & marker table below for plotting
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];
% ========================================
% define incoming files & related properties.


samplename = '#STT5442-45';
thickness = [1.85, 2.3, 4.0, 5.3];

xlim_int = [0,35]; % Hres(T)
ylim_int = [0,200]; % f(GHz)
titlename = ' $W6PN00021(t = 2.06  nm): f - H_{res}$';
xlabelname = '$f (GHz)$';
ylabelname = '$\Delta H (Oe)$';

% set tilted angle
angle_degree = 6;
angle_rad = angle_degree*3.14159265/180; % angle: degree to radian

file_format = 'In*.dat';
% ========================================

% read all data files with file_format in the current folder
files=dir(file_format);

% sort out the files in natural order
[filenames, index] = sort_nat({files.name});

%%% read files in the current folder
folder = pwd;

% open and creat the file for output
% outputloc=[folder '/' outputname];

% Open or create new file for writing. Append data to the end of the file.
% fidout=fopen(outputloc,'a+');
% ========================================
% open a figure for plotting
fig = figure();
set(fig, 'Position', [200, 100, 1000, 800]);
set(fig,'color','w');
% ========================================

i_start = 1;
i_end = numel(filenames);
h = zeros(3,1);

for i = i_start:i_end
    
fH = importdata(filenames{i});

f = fH(:,1);

Hres = fH(:,2)/10000; % Oe to T
Hres = Hres*cos(angle_rad);

% lw = fH(:,3)*cos(angle_rad);
% lw_low = fH(:,4)*cos(angle_rad);
% lw_up = fH(:,5)*cos(angle_rad);
% lw_err = (lw_up-lw_low)/2;

x = f;
y = Hres;
y_err = Hres;

% h(i) = plot(Hres,f,'color',colortable(i), 'marker',markertable(i),'markersize',20);
h(i) = errorbar(x,y,y_err,'color',colortable(i),'marker',markertable(i),'markersize',20);

hold on;

xlim(xlim_int)
ylim(ylim_int)
title(titlename,'FontSize',42,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold')

ylabel(ylabelname,'FontSize',36,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold')
xlabel(xlabelname,'FontSize',36,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold') 
% xlim([0,0.4]);
set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%g'));
set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%g'));




end
% legend(h,'1.85 nm', '2.3 nm', '4.0 nm','5.3 nm','location','northwest');
grid on;
legend('new measurement','old measurement','location','northwest');
% xmesh = linspace(0,1,100);
% fitresult = polyfit(Hres,f,1);
% ymesh = xmesh*fitresult(1)+fitresult(2);
% line(xmesh, ymesh, 'color','b','linewidth',5);
% legend(h,strsplit(num2str(thickness)));

% =========================================
% errorbar(f,lw,lw_err,'ro','markersize',20);
% xlim(xlim_int)
% ylim(ylim_int);
% % 
% % % Plotting title
% titlename = ['\DeltaH - f  of #', samplename, ', t(CoFeB) = ',.
%     num2str(samplethickness),' nm'];
% title(titlename,'FontSize',42,'FontWeight','bold');
% % 
% % % x & y labels
% xlabel(xlabelname,'FontSize',36,'FontWeight','bold');
% ylabel(ylabelname,'FontSize',36,'FontWeight','bold');
% % 
% % % axis ticks
% set(gca,'Fontsize',36,'Linewidth',3,'fontweight','bold');
% set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%g'));
% set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%g'));