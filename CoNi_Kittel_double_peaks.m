% ========================================
% Fit Kittel and Gilbert Damping 
% Extract Heff, alpha, \DeltaH0 with fixed g-factor
% ========================================
clear;
close all;
% ========================================
% Define color table & marker table below for plotting
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];



xlim_int = [0,1.0]; % Hres(T)
ylim_int = [0,4]; % f(GHz)

fieldmesh = linspace(0,1.0,100);
fmesh = linspace(0,35,100);

titlename = 'Fit Kittel (CoNi): $\frac{1}{H_{res}}(\frac{\omega}{\frac{\mu_B}{h}})^2 - H_{res}$';
xlabelname = '$H_{res} (T)$';
% ylabelname = '$\frac{1}{H_{ext}}(\frac{\omega}{\gamma_e})^2 (T)$';
ylabelname = '$\frac{1}{H_{res}}(\frac{\omega}{\frac{\mu_B}{h}})^2$';


flim_int = [0,35];
dH_int = [0, 400];
 
mu_B = 927.400*10^(-26); % J/T
h = 6.62607*10^(-34); % J*s
% ========================================

% open the following txt file to save data

%% read files in the current folder
PN = '2';


% 2.67, 2.33, 2.06, 1.86, 1.69
thickness = 0;

if strcmp(PN, '2')
    thickness = 2.67;
else 
    if strcmp(PN,'6')
        thickness = 2.33;
    else
        if strcmp(PN,'10')
            thickness = 2.06;
        else
            if strcmp(PN,'14');
                thickness = 1.86;
            else if strcmp(PN,'18');
                    thickness = 1.69;
                end
            end
        end
    end
end

fileformat = ['*_P',PN,'N1*.txt'];

files = dir(fileformat);

[filenames, index] = sort_nat({files.name});% sort out the files in natural order
files = files(index);
len_files = numel(files);
p = ones(size(len_files));
note = zeros(1,2);
% ==========================================
% create and open the file for output
outputname = 'fit_out.txt';

folder = pwd;

outputloc=[folder '/' outputname];

% Open or create new file for writing. Append data to the end of the file.
fidout=fopen(outputloc,'a+');

% fprintf(fidout,'t(nm)    Heff(T)    g    g_u dH_0    alpha    alpha_{err}\n');

fprintf(fidout,'t(nm)\tHeff1(T)\tHeff_err1(T)\tg1\tg_err1\tHeff2(T)\tHeff_err2(T)\tg2\tg_err2   \n');


%plot f-H
fig = figure();
figure(fig);
axes1 = axes('Parent',fig,'FontSize',32);
set(fig, 'Position', [80, 60, 1000, 800]);

for i = 1:len_files
    
   for k = 1:2
    
data = importdata(files(i).name);
f = data(:,1); % in GHz

% 1st peak
if k == 1
x = data(:,2)*10^(-4); % in Tesla
y = (f*10^(9)).^2./((mu_B/h)^2*x);
plotcolor = 'r';
plotmarker = 'o';

% 2nd peak
else 
x = data(:,6)*10^(-4); % in Tesla
y = (f*10^(9)).^2./((mu_B/h)^2*x);
plotcolor = 'b';
plotmarker = 's';
end

plot(x,y,[plotcolor,plotmarker],'MarkerSize',20);
title(titlename,'FontSize',42,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlim(xlim_int);
ylim auto;

ylabel(ylabelname,'FontSize',36,'FontWeight','bold', 'interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlabel(xlabelname,'FontSize',36,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
hold on;


% select data to be fitted.
[xleft,yleft]=ginput(1);
[xright,yright]=ginput(1);
 
ind=zeros(length(x),1);

    for j =1:1:length(x);
         if x(j)<xleft || x(j)>xright
            ind(j)=j;
         end               
    end
    
x_slt = x(setdiff(1:length(x),ind));
y_slt = y(setdiff(1:length(x),ind));
f_slt = f(setdiff(1:length(x),ind));

% lw_slt = lw(setdiff(1:length(x),ind));
% lw_low_slt = lw_low(setdiff(1:length(x),ind));
% lw_up_slt = lw_up(setdiff(1:length(x),ind));

% Plot out selected data
plot(x_slt,y_slt,[plotcolor,plotmarker],'MarkerSize',10);

% =========================================================
% Model with 2-parameters
testx = x_slt;
testy = y_slt;

ok_ = isfinite(testx) & isfinite(testy); %% just checking testx and testy are finite (wich are the plot we want to fit)
%exclude data exclude points from testx,testy through the logic ok_ variable. if ok_ is zero we exclude that point
x2=testx(ok_);
y2=testy(ok_);
%%the fitting domain is the whole data range
%%%fitting with a A*sqrt((x-ex)*(x-ex+xeff)) ex:exchange bias
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[-100 -100],'Upper',[100 100],...
'DiffMinChange', 1e-16,'TolFun', 1e-14 ,'MaxIter',15000,'MaxFunEvals',15000,...
'Exclude',excludedata(testx(ok_), testy(ok_),'domain', [x2(1),x2(length(x2))]));
            
len_x = length(x2);
st_A = (y2(len_x)-y2(1))/(x2(len_x)-x2(1));
st_eff = y2(len_x)/y2(len_x)-x2(len_x);
            
st_ = [st_A, st_eff];%initial condition
            
set(fo_,'Startpoint',st_);
ft_ = fittype('A*(x+xeff)',...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'A', 'xeff'});
%Fit this model using new data
[cfunP,~,~] = fit(x2,y2,ft_,fo_);
% fitting parameter --> ParamP
paramP=coeffvalues(cfunP);
%confidence of fit parameters (2 \delta region)
ciP = confint(cfunP,0.95);
% Plot this fit and write the parameters     
ymesh = paramP(1)*(fieldmesh+ paramP(2));  



line(fieldmesh, ymesh, 'color',plotcolor,'linewidth',4);
xlim(xlim_int);
ylim(ylim_int);


g=sqrt(paramP(1));
g_up = sqrt(ciP(2,1));
g_low = sqrt(ciP(1,1));
g_err = (g_up - g_low)/2;% half of error bar length
            
Heff = paramP(2);
Heff_up = ciP(2,2);
Heff_low = ciP(1,2);
Heff_err = (Heff_up-Heff_low)/2;

Kittel_equation ='$$\frac{\omega}{\gamma_e} = \sqrt{H_{res}(H_{res}+H_{eff})}$$'; 
gamma_e = '$ \gamma_e = g\frac{\mu_B}{\hbar}$';  
p1 = '$$ H_{eff} = $$';
p2=sprintf('%1.2f',Heff); % Effectrive field
p3 = '$$ \pm $$';
p4=sprintf('%1.2f T',Heff_err); % Effectrive field
text_Heff = [p1,p2,p3,p4];
             
p5 = '$$ g =  $$';
p6=sprintf('%1.2f',g); % Effectrive field
p7 = '$$ \pm $$';
p8=sprintf('%1.2f',Heff_err); % Effectrive field
text_g = [p5,p6,p7,p8];

if k == 1
    note1 = {Kittel_equation,'\#1',text_Heff,text_g};
    Heff1 = Heff;
    Heff_err1 = Heff_err;
    g1 = g;
    g_err1 = g_err;
else 
    note2 = {'\#2',text_Heff,text_g}; 
    Heff2 = Heff;
    Heff_err2 = Heff_err;
    g2 = g;
    g_err2 = g_err;
end

end
   % 
annotation(fig,'textbox',...
[0.15 0.55 0.5 0.3],...
'string',[note1,note2],'FitBoxToText','on',...
'LineStyle','none','FontSize',32,  'interpreter','latex',...
'fontsize',32,'FontWeight','bold');
grid on;


% save data
fprintf(fidout,'%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t\n',thickness,Heff1,Heff_err1, g1,g_err1,Heff2,Heff_err2, g2,g_err2);


% Save figure
    fig.PaperPositionMode = 'auto';% set image size as auto
    saveas(fig,strtok(char(filenames{i}),'.'),'png');

%    
end
close(fig);

% h_lgd = legend('1^{st} Peak','2^{nd} Peak','location','northwest');
% set(h_lgd,'fontsize',32);

% lgd_array = strread(num2str(thickness,3),'%s');
% h_lgd = legend(p,cell2mat(lgd_array),'location','southeast');
% h_lgd = legend(p,'2.67','2.33', '2.19', '2.06', '1.96', '1.86', '1.78', '1.69','location','southeast');

% saveas(fig2,'fH.png');
% fclose(fidout);

% close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



