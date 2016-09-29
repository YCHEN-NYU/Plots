% Fit Kittel & Gilber Damping for Perpendicular FMR 
% Yiyi, 09/22/2016

clear;
close all;

% Define color table & marker table below for plotting
% colortable = ['r','b','c','k','g','m','r','b',...
%     'c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^',...
    'o','s','v','^','o','s','v','^'];

colortable = lines(8);
sizeMK = 15;

mu_B = 927.400*10^(-26); % J/T
h = 6.62607*10^(-34); % J*s
factor = 10^(9)*h/mu_B;

%***************************
%***************************
% cd '/Users/yiyi/Desktop/FINAL_v2/W6PN21_P2N1_0degree';

% enter sample name
sampleName = 'W6PN'; 
if strcmp(sampleName, 'STT54')
    N = 4;
else 
    N = 3;
end

angleDegree = 6;
angleRad = angleDegree*3.14159265/180;
gMat = []; % define an empty variable to store g-factor

yleft = 7.5;
yright = 26.5;


HLowerbound = 0;
HUpperbound = 1;

fLowerbound = 0;
fUpperbound = 30;

dHLowerbound = 0;
dHUpperbound = 200;

Hlim = [HLowerbound,HUpperbound]; % Hres(T)
flim = [fLowerbound,fUpperbound]; % f(GHz)
dHlim = [dHLowerbound, dHUpperbound];
f2lim = [0, 10];

meshPoints = 100;
Hmesh = linspace(HLowerbound,fUpperbound,meshPoints);
fmesh = linspace(fLowerbound,fUpperbound,meshPoints);

% thickness = [2.67, 2.06, 1.69,2.67, 2.06, 1.69];
thickness = [1.85, 2.3, 4.0, 5.3,1.85, 2.3, 4.0, 5.3 ];
% thickness = [1.85, 1.85, 2.3, 2.3, 4.0, 4.0, 5.3, 5.3];

titlename = ['$Kittel(' sampleName  '): \frac{1}{H_{res}}(\frac{\omega}{\frac{\mu_B}{h}})^2 - H_{res}$'];
xlabelname = '$H_{res} (T)$';
ylabelname = '$\frac{1}{H_{res}}(\frac{\omega}{\frac{\mu_B}{h}})^2$';

titlename2 = ['$Gilbert(' sampleName '): \Delta H - f$'];
xlabelname2 = '$f (GHz)$';
ylabelname2 = '$\Delta H (Oe)$';


fileFormat = ['*' sampleName '*.txt'];
files=dir(fileFormat);
[filenames, index] = sort_nat({files.name});
outputname = ['output' '.txt'];
folder = pwd;
outputloc=[folder '/' outputname];
fidout=fopen(outputloc,'a+');


% fprintf(fidout,'t(nm)    Heff(T)     Heff_err(T)    ...
% g    g_err       dH0(Oe)    dH0_err(Oe)    alpha   
% alpha_err\n');

% ========================================
% open a figure for plotting
fig1 = figure();
set(fig1, 'Position', [200, 100, 1000, 800]);
set(fig1,'color','w');

fig2 = figure();
set(fig2, 'Position', [200, 100, 1000, 800]);
set(fig2,'color','w');
% ========================================


% if strcmp(sampleName, 'STT5438-41')
%     i_start = 1;
%     i_end = 4;
%     outputname = 'output_STT5438-41.txt';
%     fignameKittel = 'Kittel_STT5438-41.png';
%     fignameGilbert = 'Gilbert_STT5438-41.png';
% 
% else 
%     if strcmp(sampleName, 'STT5442-45')
%     i_start = 5;
%     i_end = 8;
%     outputname = 'output_STT5442-45.txt';
%     fignameKittel = 'Kittel_STT5442-45.png';
%     fignameGilbert = 'Gilbert_STT5442-45.png';
%     end
% end

i_start = 1;
i_end = numel(filenames);

plot1 = zeros(1, i_end);
plot2 = plot1;
legendLine = [];
legendLine_Gilbert =[];

for i = i_start:i_end
 
rawdata = importdata(filenames{i});

f = rawdata(:,1);
Hres = rawdata(:,2)/10000; % Oe to T
Hres = Hres*cos(angleRad);

lw = rawdata(:,3)*cos(angleRad);
lw_low = rawdata(:,4)*cos(angleRad);
lw_up = rawdata(:,5)*cos(angleRad);
lw_err = (lw_up-lw_low)/2;

x = Hres;
y = (f*10^(9)).^2./((mu_B/h)^2*Hres);


figure(fig1);
if(i <= N)
plot1(i) = plot(x,y, 'linestyle', 'none','color',colortable(i,:), 'marker',markertable(i),'markersize',sizeMK);
else
    plot1(i) = plot(x,y, 'linestyle', 'none','color', colortable(i-N,:), 'marker', markertable(i-N),'markersize',sizeMK,'MarkerFaceColor',colortable(i-N,:));
end
hold on;


title(titlename,'FontSize',42,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlim(Hlim);
ylim(f2lim);


ylabel(ylabelname,'FontSize',36,'FontWeight',...
    'bold', 'interpreter','latex',...
    'fontsize',42,'FontWeight','bold');
xlabel(xlabelname,'FontSize',36,'FontWeight',...
    'bold','interpreter','latex',...
    'fontsize',42,'FontWeight','bold');

set(gca,'Fontsize',36,'Linewidth',3,'fontweight','bold');
% select fitting area
% [xleft,yleft]=ginput(1);
% [xright,yright]=ginput(1);
 
fleft = 10;
fright = 28;

ind=zeros(length(f),1);

    for j =1:1:length(f);
         if f(j)< fleft || f(j)> fright
            ind(j)=j;
         end               
    end
    
x_slt = x(setdiff(1:length(f),ind));
y_slt = y(setdiff(1:length(f),ind));
f_slt = f(setdiff(1:length(f),ind));

lw_slt = lw(setdiff(1:length(f),ind));
lw_low_slt = lw_low(setdiff(1:length(f),ind));
lw_up_slt = lw_up(setdiff(1:length(f),ind));

% Plot out selected data
if i <= N
plot(x_slt,y_slt, 'linestyle', 'none','color', colortable(i,:), 'marker', markertable(i),'MarkerSize',sizeMK);
else
    plot(x_slt,y_slt, 'linestyle', 'none','color', colortable(i-N,:), 'marker', markertable(i-N),'MarkerSize',sizeMK);
end
% =========================================================
% Model with 2-parameters
testx = x_slt;
testy = y_slt;
ok_ = isfinite(testx) & isfinite(testy); 

x2=testx(ok_);
y2=testy(ok_);
 
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 -100],'Upper',[100 100],...
'DiffMinChange', 1e-16,'TolFun', 1e-14 ,'MaxIter',...
15000,'MaxFunEvals',15000,...
'Exclude',excludedata(testx(ok_), testy(ok_),'domain',...
[x2(1),x2(length(x2))]));
            
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

% ====================================================
% Model with 1-parameters (g is fixed here)
% Fit data with y = x + c model (fixed )
% f^2/((Gamma/2pi)*Hres) = Hres + Heff
% fit selected data
% 
% [xData, yData] = prepareCurveData( x_slt, y_slt );
% 
% % Set up fittype and options and starting point
% 
% ft = fittype( 'x+c', 'independent', 'x', 'dependent', 
% 'y' ); % model
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); 
% options
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
% line(Hmesh,Hmesh+fitresult.c,'linewidth',2,'color',
% colortable(i));
% ====================================================

figure(fig1);
if i <= N
plot(fmesh, cfunP(fmesh), 'color',colortable(i,:),'LineWidth',2);
else
    plot(fmesh, cfunP(fmesh), 'color',colortable(i-N,:),'LineWidth',2);
end
xlim(xlim);
ylim(ylim);

g=sqrt(paramP(1));
g_up = sqrt(ciP(2,1));
g_low = sqrt(ciP(1,1));
g_err = (g_up - g_low)/2;% half of error bar length
            
Heff = paramP(2);
Heff_up = ciP(2,2);
Heff_low = ciP(1,2);
Heff_err = (Heff_up-Heff_low)/2;

% Kittel_equation ='$$\frac{\omega}{\gamma_e}_\perp = 
% H_{res}+H_{eff}$$'; 
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
% 'string',{Kittel_equation,text_Heff,text_g},
% 'FitBoxToText','on',...
% 'LineStyle','none','FontSize',32,  'interpreter',
% 'latex',...
% 'fontsize',32,'FontWeight','bold');
%      
clear testx testy;
 
% % =============================================
figure(fig2);
testx = f_slt;
testy = lw_slt;

%%%exclude all zero points
ok1=excludedata(testx,testy,'box',[0 testx(length(testx))...
    -1e-4 1e-4]);
ok2=excludedata(testx,testy,'box',[0 10 10e-3 1]);
ok3=excludedata(testx,testy,'box',[0 testx(length(testx)) ...
    40e-3 1]);

ok_ = isfinite(testx) & isfinite(testy);  
x1=testx(ok_); %x is Frequencies in GHz
y1=testy(ok_); %y is width in Oe
            
% 'marker',markertable(i),'MarkerSize',10)
if(i <= N)
plot2(i) = errorbar(f,lw,lw_err, 'linestyle', 'none','color', colortable(i,:), 'marker', markertable(i),'MarkerSize',sizeMK);
else 
    plot2(i) = errorbar(f,lw,lw_err, 'linestyle', 'none','color', colortable(i-N,:),'marker', markertable(i-N),'MarkerSize',sizeMK,'MarkerFaceColor',colortable(i-N,:));
end
hold on;
% plot(f_slt,lw_slt,[colortable(i) markertable(i)],'MarkerSize',5)


title(titlename2,'FontSize',42,'FontWeight',...
    'bold','interpreter','latex',...
'fontsize',36,'FontWeight','bold')

ylabel(ylabelname2,'FontSize',42,'FontWeight',...
    'bold','interpreter','latex',...
'fontsize',36,'FontWeight','bold')
xlabel(xlabelname2,'FontSize',42,'FontWeight',...
    'bold','interpreter','latex',...
'fontsize',36,'FontWeight','bold') 

set(gca,'Fontsize',36,'Linewidth',3,'fontweight','bold');
xlim(flim);
ylim(dHlim);


            
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0  0],...
    'Upper',[200 100],'DiffMinChange', 1e-18,'TolFun',...
    1e-18 ,'MaxIter',15000,...
    'MaxFunEvals',15000, 'Exclude',excludedata(testx(ok_), ...
    testy(ok_),'domain', [min(x1) max(x1)]));

st_alpha = (f_slt(length(f_slt))...
    -f_slt(1))/(lw_slt(length(lw_slt))...
    -lw_slt(1));
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

gamma=2*pi*1e9*1e4/(1.758820*1e11*0.5*g);
% 1e9 from GHz, 1e4 from Oe 

dH0_low = sqrt(ciPlw(1,1))/gamma;
dH0_up = sqrt(ciPlw(2,1))/gamma;
dH0_err = (dH0_up - dH0_low)/2;

dH0 = paramPlw(1);

%alpha-factor with lower and upper bound for 95% confidence

alpha_low = sqrt(ciPlw(1,2))/gamma;
alpha_up = sqrt(ciPlw(2,2))/gamma;

alpha = paramPlw(2)/gamma;
alpha_err = (alpha_up - alpha_low)/2;
     

figure(fig2);
% Plot this fit and write the parameters
if i <= N
plot(fmesh, cfunPlw(fmesh), ...
    'color',colortable(i,:),'LineWidth',2);
else
    plot(fmesh, cfunPlw(fmesh), ...
    'color',colortable(i-N,:),'LineWidth',2);
end
% GilbertEquation ='$$\Delta H = \Delta H_0 + 
% (\frac{2 \pi}{\gamma_e})f \alpha$$'; 
          
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

fprintf(fidout,...
'%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%5.3g\t%5.3g\t\n',...
    thickness(i),Heff,Heff_err, g,g_err, dH0,...
    dH0_err, alpha, alpha_err);

legendThickness = sprintf('%1.2f',thickness(i)); 
legendSpacer = '\hspace{8 mm}';
legendGfactor = sprintf('%1.2f',g);
legendHeff = sprintf('%1.2f',Heff);


legendLine  = [legendLine; legendThickness ...
    legendSpacer legendGfactor legendSpacer legendHeff];

alpha_disp = alpha*10^3;

legendThickness_Gilbert = sprintf('%5.2f',thickness(i)); 
legendSpacer_Gilbert = '\hspace{8 mm}';
legendAlpha_Gilbert = sprintf('%5.1f',alpha_disp);
legenddH0_Gilbert = sprintf('%5.1f',dH0);


legendLine_Gilbert  = [legendLine_Gilbert; ...
    legendThickness_Gilbert legendSpacer_Gilbert legendAlpha_Gilbert legendSpacer_Gilbert legenddH0_Gilbert];
end

% 
% % **************************************************
legendHeader ='$t(nm)\hspace{6mm} g \hspace{10mm}   H_{eff}(T)$';
KittelEquation = ...
    '$$(\frac {\omega}{\gamma_e})^2_{\parallel} = H_{res}(H_{res}+H_{eff})$$'; 

% thickness = [1.85, 2.3, 4.0, 5.3 ];


figure(fig1);
% legend(plot1, '1.69 nm', '2.06 nm','2.67 nm', ...
%     'location', 'northwest');
% legend(plot1, '2.67 nm - B05', '2.06 nm - B05','1.69 nm - B05','2.67 nm - B03', '2.06 nm - B03','1.69 nm - B03', ...
%      'location', 'northwest');
%  

if strcmp(sampleName, 'STT54')

LGD1 = legend(plot1, '1.85 nm - old', '2.3 nm - old',  '4.0 nm - old','5.3 nm - old', ...
   '1.85 nm - new','2.3 nm - new',  '4.0 nm - new', '5.3 nm - new','location', 'northwest');
else 
LGD1 = legend(plot1, '2.67 nm - old', '2.06 nm - old','1.69 nm - old','2.67 nm - new', '2.06 nm - new','1.69 nm - new', ...
     'location', 'northwest');
end
% LGD1 = legend(plot1, '1.85 nm - old', '2.3 nm - old',  '4.0 nm - old','5.3 nm - old', ...
%    '1.85 nm - new','2.3 nm - new',  '4.0 nm - new', '5.3 nm - new','location', 'southeast');
set(LGD1, 'fontsize', 24);

% annotation(fig1,'textbox',...
% [0.55 0.2 0.5 0.3],...
% 'string',{KittelEquation, '\newline', legendHeader, legendLine}...
% ,'FitBoxToText','on',...
% 'LineStyle','none','FontSize',24,  'interpreter','latex',...
% 'fontsize',24,'FontWeight','bold');



% % **************************************************
% % 1.85, 2.3, 4.0, 5.3
% % **************************************************
GilbertEquation = ...
'$$\Delta H = \Delta H_0 +  (\frac{4 \pi}{\gamma_e}) \alpha f$$'; 

legendHeader_Gilbert ='$t(nm)\hspace{2mm} \alpha \times 10^{-3} \hspace{2mm}   \Delta H_{0}(Oe)$';

figure(fig2);
% 
% annotation(fig2,'textbox',...
% [0.625 0.675 0.5 0.3],...
% 'string',{GilbertEquation,'\newline', legendHeader_Gilbert,...
% legendLine_Gilbert}...
% ,'FitBoxToText','on',...
% 'LineStyle','none','FontSize',24,  'interpreter','latex',...
% 'fontsize',24,'FontWeight','bold');

% 
% 
% annotation(fig2,'textbox',...
% [0.45 0.55 0.5 0.3],...
% 'string',{GilbertEquation}, 'FitBoxToText','on',...
% 'LineStyle','none','FontSize',24, 'interpreter','latex',...
% 'fontsize',24,'FontWeight','bold');
% **************************************************


figure(fig2);
%  
if strcmp(sampleName, 'STT54')

LGD2 = legend(plot2, '1.85 nm - old', '2.3 nm - old',  '4.0 nm - old','5.3 nm - old', ...
   '1.85 nm - new','2.3 nm - new',  '4.0 nm - new', '5.3 nm - new','location', 'northwest');
else 
LGD2 = legend(plot2, '2.67 nm - old', '2.06 nm - old','1.69 nm - old','2.67 nm - new', '2.06 nm - new','1.69 nm - new', ...
     'location', 'northwest');
end
set(LGD2, 'fontsize', 24);

fignameKittel = [sampleName '_Kittel.png'];
fig1.PaperPositionMode = 'auto';% set image size as auto
saveas(fig1,  fignameKittel);

fignameGilbert = [sampleName '_Gilbert.png'];

fig2.PaperPositionMode = 'auto';% set image size as auto
saveas(fig2,  fignameGilbert);

% close all;