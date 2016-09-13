clear;
close all;
% ========================================
% Define color table & marker table below for plotting
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];

% ========================================
% Move to destination folder
cd '/Users/yiyi/Desktop/Eason_FMR/STT5438-5445A/STT5440A';
fidout=fopen('fit_fH.txt','a+'); % Load output file 
% ========================================
% constants for unit conversion
% Check "The NIST Reference on Constants, Units, and Uncertainty" for more details 
% g = 2.135;
em=1.758820*1e11; %electron charge to mass quotient 
% A = g^2*em^2/(4*pi*1e9)^2;
% gamma=2*pi*1e9*1e4/(1.758820*1e11*0.5*g);% 1e9 from GHz, 1e4 from Oe 

% ========================================

%% open the following txt file to save data

%%% read files in the current folder
files=dir('fit_parameters.txt');
[filenames, index] = sort_nat({files.name});% sort out the files in natural order
files = files(index);
len_files = numel(files);
p = ones(size(len_files));


for i = 1:len_files
data = importdata(files(i).name);
f = data(:,1); % in GHz
Hres = data(:,2);
Hres = Hres*10^(-4); % in Tesla


% Hres2 = data(:,6);
% Hres2 = Hres2*10^(-4);


x = Hres;
y = f.^2./Hres;
% y2 = f.^2./Hres2;

% plot f-H

fig = figure();
set(fig, 'Position', [200, 100, 1000, 800])
set(gcf,'color','w');
plot(Hres,y,'ro','MarkerSize',20)

hold on;

% plot(Hres2,y2,'bs','MarkerSize',20);

% xlim([0,0.4]);
% ylim([1000,1600])
xlabel('H_{res}(T)','FontSize',32,'FontWeight','bold') 
ylabel('f^2/H_{res}(GHz^2/T)','FontSize',32,'FontWeight','bold') 
title('f^2/H_{res} - H_{res} of #STT5440A, t = 4 nm','fontsize',42,'fontweight','b');

set(gca,'Fontsize',32,'Linewidth',3,'fontweight','bold');
%set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
%set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));
set(gcf,'NextPlot','add');
axes; 
set(gca,'Visible','off'); 
% set(h,'Visible','on');
%==========================


            
% select data to be fitted.
[xleft,yleft]=ginput(1);
[xright,yright]=ginput(1);
 
ind=zeros(length(x),1);

    for j =1:1:length(x);
         if x(j)<xleft || x(j)>xright
            ind(j)=j;
         end               
    end
    
x_slt=x(setdiff(1:length(x),ind));
y_slt=y(setdiff(1:length(y),ind));



testx = x_slt;
testy = y_slt;
ok_ = isfinite(testx) & isfinite(testy); %% just checking testx and testy are finite (wich are the plot we want to fit)
            %exclude data exclude points from testx,testy through the logic ok_ variable. if ok_ is zero we exclude that point
            x2=testx(ok_);
            y2=testy(ok_);
            %%the fitting domain is the whole data range
            %%%fitting with a A*sqrt((x-ex)*(x-ex+xeff)) ex:exchange bias
            fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 -5],'Upper',[2000 5],...
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



% Plot this fit and write the parameters
            fieldmesh = linspace(0,0.4,1000);
            
            paramP=coeffvalues(cfunP);
            
            fmesh = cfunP(fieldmesh);
            
            plot(fieldmesh, fmesh, 'r','LineWidth',4);
            
            %confidence of fit parameters (2 \delta region)
            ciP = confint(cfunP,0.95);
            
            % g factor
            g=sqrt(paramP(1))/em*2*2*pi*1e9;
            g_low = sqrt(ciP(1,1))/em*2*2*pi*1e9;
            g_up = sqrt(ciP(2,1))/em*2*2*pi*1e9;
            g_u = (g_up - g_low)/2;
            % Fitting Parameters
            p1=sprintf('= %5.3f GHz/T ',sqrt(paramP(1))); % slope
            p2=sprintf('H_{eff}= %5.3f T ',paramP(2)); % Effectrive field
            p3=sprintf('g value= %5.3f',g);
            p4=sprintf('g_{err} = %5.3f',g_u);% 2\delta of g factor

            annotation(fig,'textbox',...
            [0.6 0.1 0.6 0.3],...
            'String',{['\gamma /2\pi' p1],p2,p3,p4},'FitBoxToText','off',...
            'LineStyle','none','FontSize',32);
        
           % waitforbuttonpress;
            params=[sqrt(paramP(1)),g,paramP(2)];
            Heff1 = paramP(2);
            Slope1 = sqrt(paramP(1));

end