
clear;
close all;

% Define color table & marker table below for plotting
% cc = hsv(8);
% colortable = hsv(8);
% colortable = ['r','b','k','c','g','m','r','b',...
%     'k','c','g','m','r','b','c','k','g','m'];

markertable = ['o','s','v','^','o','s','v','^',...
    'o','s','v','^','o','s','v','^'];
colortable = lines(8);
sizeMK = 15;
% ***********************************
sampleName = 'W6PN';
if strcmp(sampleName, 'STT54')
    N = 4;
else
    N = 3;
end

angleDegree = 6;
angleRad = angleDegree*3.14159265/180;
% ***********************************


titlename = ['$f - H(' sampleName  ')$'];
xlabelname = '$H_{res} (T)$';
ylabelname = '$f (GHz)$';


titlename2 = ['$Gilbert(' sampleName '): \Delta H - f$'];
xlabelname2 = '$f (GHz)$';
ylabelname2 = '$\Delta H (Oe)$';


yleft = 7.5;
yright = 26.5;


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

fig1 = figure();
set(fig1, 'Position', [200, 100, 1000, 800]);
set(fig1,'color','w');
fig1.PaperPositionMode = 'auto';% set image size as auto


fig2 = figure();
set(fig2, 'Position', [200, 100, 1000, 800]);
set(fig2,'color','w');
fig2.PaperPositionMode = 'auto';% set image size as auto

plot1 = zeros(1,2*N);
plot2 = plot1;


for i = 1:2*N
    rawdata = importdata(filenames{i});
    
    f = rawdata(:,1);
    H = rawdata(:,2)/10^4*cos(angleRad);
    lw = rawdata(:,3)*cos(angleRad);
    lw_lb = rawdata(:,4)*cos(angleRad);
    lw_ub = rawdata(:,5)*cos(angleRad);
    
%     if i > 3 
%         H = H+0.032; % calibration of Hall sensor difference
%     end
    
    lw_err = (lw_ub -lw_lb)/2;
    
%% Plot1 = f-H plot 
figure(fig1);
    if i <=N
    plot1(i) = plot(H,f,'linestyle', 'none','color',colortable(i,:), ...
    'marker',markertable(i),'markersize',sizeMK);
    else 
        plot1(i) = plot(H,f,'linestyle', 'none','color',colortable(i-N,:), ...
    'marker',markertable(i-N),'markersize',sizeMK,'MarkerFaceColor',colortable(i-N,:));
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



%% Plot1 = f-H plot 
figure(fig2);
    if i <=N
    plot2(i) = errorbar(f,lw,lw_err,'linestyle', 'none','color',colortable(i,:), ...
    'marker',markertable(i),'markersize',sizeMK);
    else 
        plot2(i) = errorbar(f,lw,lw_err, 'linestyle', 'none','color',colortable(i-N,:), ...
    'marker',markertable(i-N),'markersize',sizeMK,'MarkerFaceColor',colortable(i-N,:));
    end
    hold on;
    

title(titlename2,'FontSize',42,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlim(flim);
ylim(dHlim);


ylabel(ylabelname2,'FontSize',36,'FontWeight',...
    'bold', 'interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlabel(xlabelname2,'FontSize',36,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',36,'Linewidth',3,'fontweight','bold');


end


figure(fig1);

if strcmp(sampleName, 'STT54')

LGD1 = legend(plot1, '1.85 nm - old', '2.3 nm - old',  '4.0 nm - old','5.3 nm - old', ...
   '1.85 nm - new','2.3 nm - new',  '4.0 nm - new', '5.3 nm - new','location', 'northeast');
else 
LGD1 = legend(plot1, '2.67 nm - old', '2.06 nm - old','1.69 nm - old','2.67 nm - new', '2.06 nm - new','1.69 nm - new', ...
     'location', 'east');
end

set(LGD1, 'fontsize',20);
fignameKittel = [sampleName '_fH.png'];

saveas(fig1,fignameKittel);

figure(fig2);
if strcmp(sampleName, 'STT54')

LGD2 = legend(plot2, '1.85 nm - old', '2.3 nm - old',  '4.0 nm - old','5.3 nm - old', ...
   '1.85 nm - new','2.3 nm - new',  '4.0 nm - new', '5.3 nm - new','location', 'northwest');
else 
LGD2 = legend(plot2, '2.67 nm - old', '2.06 nm - old','1.69 nm - old','2.67 nm - new', '2.06 nm - new','1.69 nm - new', ...
     'location', 'northwest');
end
set(LGD2, 'fontsize',20);

fignameGilbert = [sampleName '_lwf.png'];
saveas(fig2,fignameGilbert);

close all;
