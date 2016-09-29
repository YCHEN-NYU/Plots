% ========================================
clear;
close all;
% ========================================
% Define color table & marker table below for plotting
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];
% ========================================
% define incoming files & related properties.
% Check "The NIST Reference on Constants, Units, and Uncertainty" for more details 

mu_B = 927.400*10^(-26); % J/T
h = 6.62607*10^(-34); % J*s
% ******************************
% ******************************
samplename = 'STT.txt';
% ******************************
% ******************************

% thickness = [1.85, 2.3, 4.0, 5.3];
thickness = 2.06;
% thickness = [1.85, 4.0, 5.3];

xlim_int = [0,2]; % Hres(T)
ylim_int = [0,30]; % f(GHz)

fieldmesh = linspace(0,2,100);
fmesh = linspace(0,35,100);

titlename = 'Fit Kittel (CoFeB): $\frac{1}{H_{res}}(\frac{\omega}{\frac{\mu_B}{h}})^2 - H_{res}$';
xlabelname = '$H_{ext} (T)$';
% ylabelname = '$\frac{1}{H_{ext}}(\frac{\omega}{\gamma_e})^2 (T)$';
ylabelname = '$\frac{1}{H_{res}}(\frac{\omega}{\frac{\mu_B}{h}})^2$';

titlename2 = 'Fit Gilbert Damping: $\Delta H - f$';
xlabelname2 = '$f (GHz)$';
ylabelname2 = '$\Delta H (Oe)$';

flim_int = [0,35];
dH_int = [0, 400];

% set tilted angle
angle_degree = 6;
angle_rad = angle_degree*3.14159265/180; % angle: degree to radian

file_format = 'W6PN*.txt';

% ========================================

% read all data files with file_format in the current folder
files=dir(file_format);

% sort out the files in natural order
[filenames, index] = sort_nat({files.name});

%%% read files in the current folder
folder = pwd;

% create and open the file for output
outputname = 'fit_out.txt';

outputloc=[folder '/' outputname];

% Open or create new file for writing. Append data to the end of the file.
fidout=fopen(outputloc,'a+');

% fprintf(fidout,'t(nm)    Heff(T)    g    g_u dH_0    alpha    alpha_{err}\n');

fprintf(fidout,'t(nm)    Heff(T)     Heff_err(T)    g    g_err       dH0(Oe)    dH0_err(Oe)    alpha    alpha_err\n');

% ========================================
% open a figure for plotting
fig1 = figure();
set(fig1, 'Position', [200, 100, 1000, 800]);
set(fig1,'color','w');

fig2 = figure();
set(fig2, 'Position', [200, 100, 1000, 800]);
set(fig2,'color','w');
% ========================================
plot1 = zeros(size(filenames));
plot2 = plot1;





outputname = ['output_' samplename '.txt'];
fignameKittel = ['Kittel_' samplename '.png'];
fignameGilbert = ['Gilbert_' samplename '.png'];




i_start = 1;
i_end = numel(filenames);

for i = i_start:i_end
 
rawdata = importdata(filenames{i});

f = rawdata(:,1);
Hres = rawdata(:,2)/10000; % Oe to T
Hres = Hres*cos(angle_rad);

lw = rawdata(:,3)*cos(angle_rad);
lw_low = rawdata(:,4)*cos(angle_rad);
lw_up = rawdata(:,5)*cos(angle_rad);
lw_err = (lw_up-lw_low)/2;

x = Hres;
y = (f*10^(9)).^2./((mu_B/h)^2*Hres);

figure(fig1);
plot1(i) = plot(x,y,'color',colortable(i), 'marker',markertable(i),'markersize',20);
hold on;

title(titlename,'FontSize',42,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlim(xlim_int);
ylim auto;

ylabel(ylabelname,'FontSize',36,'FontWeight','bold', 'interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlabel(xlabelname,'FontSize',36,'FontWeight','bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');

% =========================================================
% easy to get tricks error in the following code
% set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%g'));
% set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%g'));
% =========================================================
% =========================================================
% select fitting area
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

lw_slt = lw(setdiff(1:length(x),ind));
lw_low_slt = lw_low(setdiff(1:length(x),ind));
lw_up_slt = lw_up(setdiff(1:length(x),ind));

% Plot out selected data
plot(x_slt,y_slt,'color',colortable(i),'marker',markertable(i),'MarkerSize',10);

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
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0],'Upper',[100 100],...
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
                        
plot(fieldmesh, cfunP(fieldmesh), 'color',colortable(i),'LineWidth',4);
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
% 
% annotation(fig1,'textbox',...
% [0.15 0.55 0.5 0.3],...
% 'string',{Kittel_equation,text_Heff,text_g},'FitBoxToText','on',...
% 'LineStyle','none','FontSize',32,  'interpreter','latex',...
% 'fontsize',32,'FontWeight','bold');
%      
clear testx testy;
                   
% =========================================================
% Model with 1-parameters (g is fixed here)
% Fit data with y = x + c model (fixed )
% f^2/((Gamma/2pi)*Hres) = Hres + Heff
% fit selected data
% 
% [xData, yData] = prepareCurveData( x_slt, y_slt );
% 
% % Set up fittype and options and starting point
% 
% ft = fittype( 'x+c', 'independent', 'x', 'dependent', 'y' ); % model
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); % options
% opts.Display = 'Off'; 
% 
% opts.StartPoint = yleft - xleft; % starting point
% 
% [fitresult, ~] = fit( xData, yData, ft, opts );
% 
% ciP = confint(fitresult,0.95);
% 
% Heff = fitresult.c; % in Tesla
% 
% Heff_err = (ciP(2,1)-ciP(1,1))/2;
% 
% % plot fitting results
% line(fieldmesh,fieldmesh+fitresult.c,'linewidth',2,'color',colortable(i));


% =========================================================

% 
% 
% % =======================================================================
figure(fig2);

testx = f_slt;
testy = lw_slt;

%%%exclude all zero points
ok1=excludedata(testx,testy,'box',[0 testx(length(testx)) -1e-4 1e-4]);
ok2=excludedata(testx,testy,'box',[0 10 10e-3 1]);
ok3=excludedata(testx,testy,'box',[0 testx(length(testx)) 40e-3 1]);

ok_ = isfinite(testx) & isfinite(testy); %&ok1&ok2&ok3; %% just checking testx and testy are finite (wich are the plot we want to fit)
%exclude data exclude points from testx,testy through the logic ok_ variable. if ok_ is zero we exclude that point
x1=testx(ok_); %x is Frequencies in GHz
y1=testy(ok_); %y is width in Oe
            
% plot(testx(ok_),testy(ok_),'color',colortable(i),'marker',markertable(i),'MarkerSize',10)

plot2(i) = errorbar(f,lw,lw_err,'color',colortable(i),'marker',markertable(i),'MarkerSize',20);

hold on;
plot(f_slt,lw_slt,'color',colortable(i),'marker',markertable(i),'MarkerSize',10)


title(titlename2,'FontSize',42,'FontWeight','bold','interpreter','latex',...
'fontsize',32,'FontWeight','bold')

ylabel(ylabelname2,'FontSize',36,'FontWeight','bold','interpreter','latex',...
'fontsize',32,'FontWeight','bold')
xlabel(xlabelname2,'FontSize',36,'FontWeight','bold','interpreter','latex',...
'fontsize',32,'FontWeight','bold') 

set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');

% set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%g'));
% set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%g'));

xlim(flim_int);
ylim(dH_int);


            
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0  0 ],...
    'Upper',[200 100],'DiffMinChange', 1e-18,'TolFun', 1e-16 ,'MaxIter',15000,...
    'MaxFunEvals',15000, 'Exclude',excludedata(testx(ok_), testy(ok_),'domain', [min(x1) max(x1)]));

st_alpha = (f_slt(length(f_slt))-f_slt(1))/(lw_slt(length(lw_slt))-lw_slt(1));
st_lw = lw_slt(length(lw_slt))-lw_slt(1);


st_ = [100, 10];

set(fo_,'Startpoint',st_);
            
%alpha in Tesla/GHz
ft_ = fittype('A0+alpha*x',...
'dependent',{'y'},'independent',{'x'},...
 'coefficients',{'A0', 'alpha'});

%  Fit this model using new data
[cfunPlw,gof,output] = fit(testx(ok_),testy(ok_),ft_,fo_);

paramPlw=coeffvalues(cfunPlw);
            
%confidence of fit parameters (2 \delta region)
ciPlw = confint(cfunPlw,0.95);

gamma=2*pi*1e9*1e4/(1.758820*1e11*0.5*g);% 1e9 from GHz, 1e4 from Oe 

dH0_low = sqrt(ciPlw(1,1))/gamma;
dH0_up = sqrt(ciPlw(2,1))/gamma;
dH0_err = (dH0_up - dH0_low)/2;

dH0 = paramPlw(1);

%alpha-factor with lower and upper bound for 95% confidence

alpha_low = sqrt(ciPlw(1,2))/gamma;
alpha_up = sqrt(ciPlw(2,2))/gamma;

alpha = paramPlw(2)/gamma;
alpha_err = (alpha_up - alpha_low)/2;
           
% Plot this fit and write the parameters
plot(fmesh, cfunPlw(fmesh), 'color',colortable(i),'LineWidth',4);


Gilbert_linewidth ='$$\Delta H = \Delta H_0 +  (\frac{2 \pi}{\gamma_e})f \alpha$$'; 
          
p9 = '$$ \Delta H = $$';
p10=sprintf('%1.3f',dH0); % Effectrive field
p11 = '$$ \pm $$';
p12=sprintf('%1.3f Oe',dH0_err); % Effectrive field

text_dH0 = [p9,p10,p11,p12];
             
p13 = '$$ \alpha =  $$';
p14=sprintf('%1.5f',alpha); % Effectrive field
p15 = '$$ \pm $$';
p16 =sprintf('%1.5f',alpha_err); % Effectrive field
text_alpha = [p13,p14,p15,p16];

% annotation(fig2,'textbox',...
% [0.15 0.55 0.5 0.3],...
% 'string',{text_dH0,text_alpha}, 'FitBoxToText','on',...
% 'LineStyle','none','FontSize',32, 'interpreter','latex',...
% 'fontsize',32,'FontWeight','bold');


fprintf(fidout,'%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%5.3g\t%5.3g\t\n',thickness(i),Heff,Heff_err, g,g_err, dH0, dH0_err, alpha, alpha_err);

end

% 
% annotation(fig1,'textbox',...
% [0.15 0.55 0.5 0.3],...
% 'string',{'CoFeB (1.85 nm)',text_Heff,text_g},'FitBoxToText','on',...
% 'LineStyle','none','FontSize',32,  'interpreter','latex',...
% 'fontsize',32,'FontWeight','bold');
% 
% annotation(fig2,'textbox',...
% [0.3 0.55 0.5 0.3],...
% 'string',{text_dH0,text_alpha},'FitBoxToText','on',...
% 'LineStyle','none','FontSize',32,  'interpreter','latex',...
% 'fontsize',32,'FontWeight','bold');
% 
% fig1.PaperPositionMode = 'auto';% set image size as auto
% saveas(fig1,  fignameKittel);
% 
% fig2.PaperPositionMode = 'auto';% set image size as auto
% saveas(fig2,  fignameGilbert);
% 
% % clear all;
% close all;



% 
% 
% 
% annotation(fig1,'textbox',...
% [0.5 0.2 0.5 0.3],...
% 'string',{Kittel_equation,gamma_e, },'FitBoxToText','on',...
% 'LineStyle','none','FontSize',32,  'interpreter','latex',...
% 'fontsize',32,'FontWeight','bold');

% Gilbert_equation ='$\Delta H = \Delta H_0 +  (\frac{2 \pi}{\gamma_e})f \alpha$'; 
% 


% legend(plot1,'+5','+10','-5','-10', 'location','northwest');
% legend(plot2,'+5','+10','-5','-10', 'location','northwest');

% legend(plot1,'2.3 nm','4.0 nm','5.3 nm','location','northwest');
% legend(plot2,'2.3 nm','4.0 nm','5.3 nm','location','northwest');

% legend(plot1,'1.85 nm','2.3 nm','4.0 nm','5.3 nm','location','northwest');
% legend(plot2,'1.85 nm', '2.3 nm','4.0 nm','5.3 nm','location','northwest');
