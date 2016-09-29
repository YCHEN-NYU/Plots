markertable = ['o','s','v','^','o','s','v','^',...
    'o','s','v','^','o','s','v','^'];
colortable = lines(8);
sizeMK = 15;
titlename = '$H - V$';
xlabelname = '$V (Volt) $';
ylabelname = '$H (mT)$';

fileFormat = 'calH.lvm';
files=dir(fileFormat);
[filenames, index] = sort_nat({files.name});

fig = figure();
set(fig, 'Position', [200, 100, 1000, 800]);
set(fig,'color','w');
fig.PaperPositionMode = 'auto';% set image size as auto
Vside = zeros(1401,1);
Hside = Vside;
Vcenter = Vside;
Hcenter = Vside;
for i = 1:numel(filenames)
    rawdata = importdata(filenames{i});
    
    V = rawdata(2:end,2);
%     V = V;
    H2000 = rawdata(2:end,3)*1000;
    H450 = rawdata(2:end,4)*1000;
%     if i < 3
%         Vcenter = Vcenter + V;
%         Hcenter = Hcenter + H;
%     else
%         Vside = Vside + V;
%         Hside = Hside + H;

% 
plot1 = line(V, H2000,'color', colortable(i,:),'linewidth', 5);
plot2 = line(V, H450,'color', colortable(i+1,:),'linewidth', 5);
    end
%      plot3 = line(V, (H450 - H2000)*1000,'color', colortable(i,:),'linewidth', 5);

        
    hold on;    
    
grid on;
% end

Vavg = (Vside + Vcenter)/4;
dH = Hside/2 - Hcenter/2;
dH = dH*1000; % Telsa to mTesla
% line(V, H2000,'color', colortable(i,:),'linewidth', 2);
% line(V, H450,'color', colortable(i,:),'linewidth', 2);

title(titlename,'FontSize',42,'FontWeight',...
'bold','interpreter','latex',...
'fontsize',42,'FontWeight','bold');

%xlim([0,1.2]);

ylabel(ylabelname,'FontSize',36,'FontWeight',...
'bold', 'interpreter','latex',...
'fontsize',42,'FontWeight','bold');
xlabel(xlabelname,'FontSize',36,'FontWeight',...
'bold','interpreter','latex',...
'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',36,'Linewidth',3,'fontweight','bold');
box on;

% lgd = legend('CenterLeft', 'CenterRight', 'SideLeft', 'SideRight', 'location', 'northwest');

lgd = legend([plot1, plot2], 'Keithley 2000','LS-450', 'location','northwest');
% lgd = legend(plot3, '$ \Delta H = LS450-Keithley2000$', 'location','northwest');

 set(lgd, 'FontSize',36,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
