
% Define color table & marker table below for plotting
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];


%%% read files in the current folder
folder = pwd;
% read all data file ended with *MHz.dat in the current folder
fileformat = '*W6*.txt';
files=dir(fileformat);
% sort out the files in natural order
[filenames, ~] = sort_nat({files.name});
% open and creat the file for output
outputloc=[folder '/' outputname];
% Open or create new file for writing. Append data to the end of the file.
fidout=fopen(outputloc,'a+');

% thickness = [1.85, 2.3, 4.0, 5.3];




fig = figure();
set(fig, 'Position', [200, 100, 1000, 800]);
set(fig,'color','w');

xlim_int = [0,2]; % Hres(T)
ylim_int = [0,25]; % f(GHz)

titlename = 'Perpendicular FMR of STT5442-45';
xlabelname = '$H (T)$';
ylabelname = '$f (GHz)$';



i_start = 1;
i_end = numel(filenames);

for i= 5:8;
mat = importdata(filenames{i});
x = mat(:,8)*1000;
y = mat(:,1);
h = plot(x,y,[colortable(i),markertable(i)],'markersize',20);
title(titlename,'FontSize',42,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
% xlim(xlim_int);
% ylim(ylim_int);

ylabel(ylabelname,'FontSize',36,'FontWeight','bold', 'interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlabel(xlabelname,'FontSize',36,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
hold on;
end
% legend('STT5438(CoFeB): 1.85 nm', 'STT5439(CoFeB): 2.3 nm','STT5440(CoFeB): 4.0 nm','STT5441(CoFeB): 5.3 nm',...
%     'STT5442(CoFeB): 1.85 nm', 'STT5443(CoFeB): 2.3 nm','STT5444(CoFeB): 4.0 nm','STT5445(CoFeB): 5.3 nm', 'location','northwest');
    legend('STT5442(CoFeB): 1.85 nm', 'STT5443(CoFeB): 2.3 nm','STT5444(CoFeB): 4.0 nm','STT5445(CoFeB): 5.3 nm', 'location','northwest');
%     legend('STT5438(CoFeB): 1.85 nm', 'STT5439(CoFeB): 2.3 nm','STT5440(CoFeB): 4.0 nm','STT5441(CoFeB): 5.3 nm', 'location','northwest');

% thickness = [1.85, 2.3, 4.0, 5.3];
