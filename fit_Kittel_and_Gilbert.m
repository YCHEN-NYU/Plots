% ========================================
% Fit Kittel and Gilbert Damping 
% Extract Heff, alpha, \DeltaH0 with fixed g-factor
% Extract g-factor separately
% ========================================
clear;
close all;
% ========================================
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];

% define constants for unit conversion
em=1.758820*1e11;
g = 2.135;
A = g^2*em^2/(4*pi*1e9)^2;
gamma=2*pi*1e9*1e4/(1.758*1e11*0.5*g);%1e9 from GHz, 1e4 from Oe 
% ========================================
% thickness matrix(nm) 
t15 =[2.67, 2.33, 2.06,1.86,1.69];
t16 = [2.67, 2.06, 1.69];
t17 = [2.67, 2.06, 1.69];
t18 =[2.67, 2.33, 2.06 1.86 1.69];
t19 = [2.67, 2.06, 1.69];
t20 = [2.67, 2.06, 1.69];
t21 = [2.67,2.33, 2.19, 2.06, 1.96, 1.86, 1.78, 1.69]';

% inhomogeneous broadening interval
hint = [50,400]; % #15 - #20
href = [0,150];% #21
% 
% cd '/Users/yiyi/Desktop/f-H/#18/18_P2N1';
thickness = t21; % thickness
h = href;% # \DeltaH intervals

fidout=fopen('21_fit_fH.txt','a+');
fieldmesh = linspace(0,0.8,100);

% ========================================

%% open the following txt file to save data

%%% read files in the current folder
files=dir('*N1.txt');
[filenames, index] = sort_nat({files.name});% sort out the files in natural order
files = files(index);
len_files = numel(files);
p = ones(size(len_files));
fig2 = figure();

set(fig2, 'Position', [80, 60, 1000, 800])

for i = 1:len_files
data = importdata(files(i).name);
f = data(:,1);
Hres = data(:,2);
Hres = Hres*10^(-4);
% lineswidth;
lw = data(:,3);
lw_low = data(:,4);
lw_up = data(:,5);
lw_err = (lw_up-lw_low)/2;
x = Hres;

y =f.^2./(A*Hres);
% plot f-H

fig1 = figure();
figure(fig1);
axes1 = axes('Parent',fig1,'FontSize',32);

set(fig1, 'Position', [80, 60, 1200, 500])

subplot(1,2,1);
plot(x,y,'color',colortable(i),'marker',markertable(i),'MarkerSize',20)

% xlim ([0,0.8]);
% ylim ([0.5,1.8]);

xlabel('H_{res} (T)','FontSize',24)
ylabel('f^2/(A*H_{res}) (T), A = (\gamma\mu_0/(2\pi))^2','FontSize',24);
set(gca,'fontsize',24);
hold on

% select data to be fitted.
[Hleft,fleft]=ginput(1);
[Hright,fright]=ginput(1);
 
ind=zeros(length(x),1);
j = 0;
    for j =1:1:length(x);
         if x(j)<Hleft || x(j)>Hright
            ind(j)=j;
         end               
    end
    
x_slt=x(setdiff(1:length(x),ind));
y_slt=y(setdiff(1:length(y),ind));
plot(x_slt,y_slt,'color',colortable(i),'marker','.','MarkerSize',10)

f_slt = f(setdiff(1:length(f),ind));
lw_slt = lw(setdiff(1:length(lw),ind));
lw_err_slt= lw_err(setdiff(1:length(lw_err),ind));
% lw_up_slt= lw_up(setdiff(1:length(lw_up),ind));


%%

[xData, yData] = prepareCurveData( x_slt, y_slt );

% Set up fittype and options.
ft = fittype( 'x+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

opts.StartPoint = Hleft - fleft;

[fitresult, gof] = fit( xData, yData, ft, opts );
line(fieldmesh,fieldmesh+fitresult.c,'linewidth',2,'color',colortable(i));
Heff = fitresult.c; % in Tesla

ciP = confint(fitresult,0.95);
Heff_err = (ciP(2,1)-ciP(1,1))/2;



%%
testx = x_slt;
testy = y_slt*A;
ok_ = isfinite(testx) & isfinite(testy); %% just checking testx and testy are finite (wich are the plot we want to fit)
            %exclude data exclude points from testx,testy through the logic ok_ variable. if ok_ is zero we exclude that point
            x2=testx(ok_);
            y2=testy(ok_);
            %%the fitting domain is the whole data range
            %%%fitting with a A*sqrt((x-ex)*(x-ex+xeff)) ex:exchange bias
            fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 -5.00],'Upper',[2000 5.0],...
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
%             Fit this model using new data
            [cfunP,gof,output] = fit(x2,y2,ft_,fo_);
          
%             fit parameters
            paramP2=coeffvalues(cfunP);
            
%             confidence of fit parameters (2 \delta region)
            ciP = predint(cfunP,0.95);
            
%             g factor with lower and upper bound for 95% confidence
            gfactor=sqrt(paramP2(1))*2*2*pi*1e9/(em);
            gfactor_low = sqrt(ciP(1))*2*2*pi*1e9/(em);
            gfactor_up = sqrt(ciP(2))*2*2*pi*1e9/(em);
            gfactor_err = (gfactor_up - gfactor_low)/2;
%      


%%
subplot(1,2,2);


errorbar(f,lw,lw_err,'color',colortable(i),'marker',markertable(i),'MarkerSize',20)
xlabel('f(GHz)','FontSize',24)
ylabel('\DeltaH (Oe)','FontSize',24)
set(gca,'fontsize',24);
hold on;
freqmesh = linspace(0,30,100);
% ============================
% select region in \DeltaH - f curve to be fitted.
fprintf('click on the LEFT!');
[f1,dH1]=ginput(1);

fprintf('click on the RIGHT!')
[f2,dH2]=ginput(1);
 
ind=zeros(length(f),1);
j = 0;
    for j =1:1:length(f);
         if f(j)<f1 || f(j)>f2
            ind(j)=j;
         end               
    end
    
f_chosen=f(setdiff(1:length(f),ind));
lw_chosen=lw(setdiff(1:length(lw),ind));









% =============================





% fit with selected region from Kittel fitting



testx = f_chosen;
testy = lw_chosen;

%%%exclude all zero points
           ok1=excludedata(testx,testy,'box',[0 testx(length(testx)) -1e-4 1e-4]);
           ok2=excludedata(testx,testy,'box',[0 10 10e-3 1]);
           ok3=excludedata(testx,testy,'box',[0 testx(length(testx))
           40e-3 1]);

            ok_ = isfinite(testx) & isfinite(testy); %&ok1&ok2&ok3; %% just checking testx and testy are finite (wich are the plot we want to fit)
            %exclude data exclude points from testx,testy through the logic ok_ variable. if ok_ is zero we exclude that point
            x1=testx(ok_); %x is Frequencies in GHz
            y1=testy(ok_); %y is width in Oe
            
% plot(testx(ok_),testy(ok_),'color',colortable(i),'marker',markertable(i),'MarkerSize',10)
plot(x1,y1,'color',colortable(i),'marker',markertable(i),'MarkerSize',10)



xlim([0,30]);
ylim(h);


            
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[-100  -100 ],...
    'Upper',[200 100],'DiffMinChange', 1e-16,'TolFun', 1e-14 ,'MaxIter',15000,...
    'MaxFunEvals',15000, 'Exclude',excludedata(testx(ok_), testy(ok_),'domain', [min(x1) max(x1)]));
st_ = [0, 0];

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
            
%alpha-factor with lower and upper bound for 95% confidence

alpha_low = sqrt(ciPlw(1,2))/gamma;
alpha_up = sqrt(ciPlw(2,2))/gamma;

alpha = paramPlw(2)/gamma;
alpha_err = (alpha_up - alpha_low)/2;
           
% Plot this fit and write the parameters
plot(freqmesh, cfunPlw(freqmesh), 'color',colortable(i),'LineWidth',2);

deltaH0 = paramPlw(1);  
p1=sprintf('= %5.3g Oe',paramPlw(1));
p2=sprintf('= %5.3g Oe/GHz',paramPlw(2));
p3=sprintf('= %5.3g ',alpha);
p4 = sprintf('= %5.3g',alpha_err);

annotation(fig1,'textbox',...
[0.7 0.2 0.4 0.2],...
'String',{['slope' p2],['\DeltaH_0' p1],['\alpha' p3],['\alpha_u' p4] },'FitBoxToText','off',...
'LineStyle','none','FontSize',12);
        

% Save fitting parameters
fprintf(fidout,'%.3g\t%.3g\t%.3g\t%3g\t%3g\t%.5g\t%.5g\t%.3g\n',thickness(i),gfactor,gfactor_err,Heff,Heff_err,alpha,alpha_err,deltaH0);

saveas(fig1,strtok(char(filenames{i}),'.'),'png')
close(fig1);


%% additional plot of f-H

figure(fig2);

hold on;
p(i) = plot(x,sqrt(y.*x*A),'color',colortable(i),'marker',markertable(i),'markersize',5);
hold on;

line(fieldmesh,sqrt(A*fieldmesh.*(fieldmesh+fitresult.c)),'linewidth',2,'color',colortable(i));
xlim([0,0.8]);
ylim([0,30]);
xlabel('H_{res} (T)','FontSize',24)
ylabel('f(GHz)','FontSize',24)
set(gca,'fontsize',24);

end

% lgd_array = strread(num2str(thickness,3),'%s');
% h_lgd = legend(p,cell2mat(lgd_array),'location','southeast');
h_lgd = legend(p,'2.67','2.33', '2.19', '2.06', '1.96', '1.86', '1.78', '1.69','location','southeast');

set(h_lgd,'fontsize',20);
saveas(fig2,'fH.png');
fclose(fidout);

close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



