%% plot out Heff-t, alpha-t, deltaalpha-t, deltaH0 - t for CoFeB

clear;
close all;
% 15
M15 = ...
[2.67	2.14	0.0257	1.19771	0.00331835	0.015241	0.00047518	71.9
2.33	2.18	0.0154	1.13371	0.00333205	0.01806	0.00029864	67
2.06	2.16	0.0741	1.10191	0.00690643	0.020714	0.0013189	85.8
1.86	2.38	0.0394	1.01091	0.0146707	0.029275	0.0012163	77.3
1.69	2.35	0.0536	0.967458	0.00961838	0.035738	0.0010226	70.1];


% 16 
M16 = [
2.67	2.16	0.00761	1.22378	0.0023454	0.015243	0.00046002	78.9
2.06	2.2	0.0157	1.13768	0.00341533	0.021693	0.00027941	85.4
1.69	2.21	0.0272	1.04064	0.00598652	0.023999	0.00062823	130];

% 17 
M17=[2.67	2.25	0.0146	1.21419	0.00781605	0.015453	0.00025314	86.5
2.06	2.21	0.018	1.1403	0.00651777	0.019863	0.00060642	125
1.69	2.45	0.0648	1.0006	0.0174656	0.026987	0.00081135	133];


% 18 
M18 =[2.67	2.18	0.0108	1.12916	0.00317728	0.0097916	0.00030988	57.9
2.33	2.28	0.0176	1.0556	0.0138511	0.011882	0.0045221	54.1
2.06	2.21	0.0149	1.00866	0.00617138	0.011976	0.00062077	82.5
1.86	2.38	0.0509	0.928548	0.0229622	0.013065	0.00099221	81.3
1.69	2.21	0.0141	0.881871	0.00683383	0.013162	0.00063582	99];


% 19 
M19 =[2.67	2.14	0.0055	1.17454	0.000828908	0.010276	0.00068422	51.6
2.06	2.14	0.0115	1.04885	0.00244604	0.0097062	0.00041663	104];
% 1.69	2.32	0.0993	0.924038	0.00977874	0.045094	0.00011566	-100];

% 20
M20 =[
2.67	2.19	0.0136	1.19602	0.00474242	0.0085496	0.00050523	92.7];
% 2.06	2.37	0.0353	1.0252	0.0126782	0.0053029	0.0012552	139
% 1.69	2.33	0.0352	0.877494	0.0106709	-0.0022514	0.001557	200];

% 21 
M21 =[
2.67	2.13	0.00616	0.931519	0.00129031	0.007709	0.0003891	30.4
2.33	2.13	0.00171	0.834422	0.00069942	0.0074563	0.00011804	20.7
2.19	2.13	0.00308	0.783515	0.000926891	0.0078877	0.00018865	23.9
2.06	2.1	0.00467	0.742056	0.00249796	0.0083932	0.00066826	31.1
1.96	2.15	0.00235	0.682321	0.00119793	0.0086561	7.99E-05	21.4
1.86	2.15	0.00324	0.647604	0.00113804	0.0081803	0.00014222	31.9
1.78	2.14	0.00324	0.597657	0.000699396	0.0090315	0.00017357	26.3
1.69	2.12	0.00555	0.552821	0.00145469	0.011579	0.00027405	25.2];

% 21 
M21_ref =[
    2.67	2.13	0.00616	0.931519	0.00129031	0.007709	0.0003891	30.4
2.33	2.13	0.00171	0.834422	0.00069942	0.0074563	0.00011804	20.7
2.06	2.1	0.00467	0.742056	0.00249796	0.0083932	0.00066826	31.1
1.86	2.15	0.00324	0.647604	0.00113804	0.0081803	0.00014222	31.9
1.69	2.12	0.00555	0.552821	0.00145469	0.011579	0.00027405	25.2];




color_table = ['r','b','c','k','g','r','b'];
marker_table = ['o','s','v','+','.','*','^'];

%%
% ==========================================================

% Heff - t
% 1    2    3    4    5        6        7            8
% t    g    gerr Heff Heff_err alpha    alpha_err    deltaH0

% mult = 1;
% a = 1;
% b =4;
% %  
% N = 7;
% h = zeros(1,N);
% 
% fig1 = figure();
% fig1.PaperPositionMode = 'auto';% set image size as auto
% 
% h(1) = plot(M15(:,a),M15(:,b)*mult,'Marker','v','MarkerSize',10,'color','b');
% hold on;
% set(fig1, 'Position', [200, 100, 800, 600])
% 
% h(2) = plot(M16(:,a),M16(:,b)*mult,'Marker','+','MarkerSize',10,'color','c');
% h(3) = plot(M17(:,a),M17(:,b)*mult,'Marker','.','MarkerSize',10,'color','k');
% h(4) = plot(M18(:,a),M18(:,b)*mult,'Marker','s','MarkerSize',10,'color','m');
% h(5) =  plot(M19(:,a),M19(:,b)*mult,'Marker','^','MarkerSize',10,'color','g');
% h(6) =  plot(M20(:,a),M20(:,b)*mult,'Marker','*','MarkerSize',10,'color','r');
% h(7) =  plot(M21(:,a),M21(:,b)*mult,'Marker','o','MarkerSize',10,'color','r');
% 
%  xlabel('CoFeB t(nm)','FontSize',32);
%  ylabel('H_{eff}(T)','FontSize',32);
%  title('CoFeB H_{eff} - t','FontSize',36);
% %  set(gca,'fontsize',32);
%  set(gca,'Fontsize',32,'Linewidth',3,'fontweight','bold');
% % set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.2f'),'fontweight','bold');
% set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%1.1f'),'fontweight','bold');
% ylim([0.39,1.41]);
% h_lgd = legend(h,'15','16', '17', '18', '19', '20', '21','location','east');

% ==========================================================
% alpha-t
% 1    2    3    4    5        6        7            8
% t    g    gerr Heff Heff_err alpha    alpha_err    deltaH0
% 
% mult = 1000;
% a = 1;
% b =6;
% c = 7;  
% N = 7;
% h = zeros(1,N);
% 
% fig1 = figure();
% fig1.PaperPositionMode = 'auto';% set image size as auto
% 
% h(1) = plot(M15(:,a),M15(:,b)*mult,'Marker','v','MarkerSize',10,'color','b');
% hold on;
% set(fig1, 'Position', [200, 100, 800, 600])
% 
% %  errorbar(M21_0420(:,a),M21_0420(:,b)*mult,M21_0420(:,c)*mult,'Marker','s','MarkerSize',10,'color','b');
% 
% h(1) = errorbar(M15(:,a),M15(:,b)*mult,M15(:,c)*mult,'Marker','v','MarkerSize',10,'color','b');
% h(2) = errorbar(M16(:,a),M16(:,b)*mult,M16(:,c)*mult,'Marker','+','MarkerSize',10,'color','c');
% h(3) = errorbar(M17(:,a),M17(:,b)*mult,M17(:,c)*mult,'Marker','.','MarkerSize',10,'color','k');
% h(4) = errorbar(M18(:,a),M18(:,b)*mult,M18(:,c)*mult,'Marker','s','MarkerSize',10,'color','m');
% h(5) = errorbar(M19(:,a),M19(:,b)*mult,M19(:,c)*mult,'Marker','^','MarkerSize',10,'color','g');
% h(6) = errorbar(M20(:,a),M20(:,b)*mult,M20(:,c)*mult,'Marker','*','MarkerSize',10,'color','r');
% h(7) = errorbar(M21(:,a),M21(:,b)*mult,M21(:,c)*mult,'Marker','o','MarkerSize',10,'color','r');
% 
% xlabel('CoFeB t(nm)','FontSize',32);
% ylabel('\alpha x 10^{-3}','FontSize',32);
% 
% set(gca,'Fontsize',32,'Linewidth',3,'fontweight','bold');
% % set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.2f'),'fontweight','bold');
% % set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%1.1f'),'fontweight','bold');
% % ylim([0.39,1.41]);
% h_lgd = legend(h,'15','16', '17', '18', '19', '20', '21','location','east');
%%
% ==========================================================
% alpha-t
% 1    2    3    4    5        6        7            8
% t    g    gerr Heff Heff_err alpha    alpha_err    deltaH0
% 
% mult = 1000;
% a = 1;
% b =6;
% c = 7; 
% 
% % N = 7;
% % h = zeros(1,N);
% % 
% fig1 = figure();
% fig1.PaperPositionMode = 'auto';% set image size as auto
% set(fig1, 'Position', [200, 100, 800, 600])

% 
% h(1) = plot(M15(:,a),(M15(:,b)+M21_ref(:,b))*mult,'Marker','v','MarkerSize',10,'color','b');
% h(2) = plot(M18(:,a),(M18(:,b)+M21_ref(:,b))*mult,'Marker','v','MarkerSize',10,'color','b');

% hold on;
% 
% 
% h(1) = errorbar(M15(:,a),(M15(:,b)-M21_ref(:,b))*mult,(M15(:,c)+M21_ref(:,c))*mult,'Marker','v','MarkerSize',10,'color','b');
% hold on;
% h(2) = errorbar(M18(:,a),(M18(:,b)-M21_ref(:,b))*mult,(M18(:,c)+M21_ref(:,c))*mult,'Marker','+','MarkerSize',10,'color','c');
% 
% xlabel('CoFeB t(nm)','FontSize',32);
% ylabel('\Delta\alpha x 10^{-3}','FontSize',32);
% % 
% set(gca,'Fontsize',32,'Linewidth',3,'fontweight','bold');
% % % set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.2f'),'fontweight','bold');
% % % set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%1.1f'),'fontweight','bold');
% % % ylim([0.39,1.41]);
% h_lgd = legend(h,'15-21','18-21','location','east');


%%
% ==========================================================

% Heff - t
% 1    2    3    4    5        6        7            8
% t    g    gerr Heff Heff_err alpha    alpha_err    deltaH0

mult = 1;
a = 1;
b = 8;
%  
N = 7;
h = zeros(1,N);

fig1 = figure();
fig1.PaperPositionMode = 'auto';% set image size as auto
set(fig1, 'Position', [200, 100, 800, 600])

h(1) = plot(M15(:,a),M15(:,b)*mult,'Marker','v','MarkerSize',10,'color','b');
hold on;

h(2) = plot(M16(:,a),M16(:,b)*mult,'Marker','+','MarkerSize',10,'color','c');
h(3) = plot(M17(:,a),M17(:,b)*mult,'Marker','.','MarkerSize',10,'color','k');
h(4) = plot(M18(:,a),M18(:,b)*mult,'Marker','s','MarkerSize',10,'color','m');
h(5) =  plot(M19(:,a),M19(:,b)*mult,'Marker','^','MarkerSize',10,'color','g');
h(6) =  plot(M20(:,a),M20(:,b)*mult,'Marker','*','MarkerSize',10,'color','r');
h(7) =  plot(M21(:,a),M21(:,b)*mult,'Marker','o','MarkerSize',10,'color','r');

 xlabel('CoFeB t(nm)','FontSize',32);
 ylabel('\DeltaH_0(Oe)','FontSize',32);
 title('CoFeB \DeltaH_0 - t','FontSize',36);
%  set(gca,'fontsize',32);
 set(gca,'Fontsize',32,'Linewidth',3,'fontweight','bold');
% set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.2f'),'fontweight','bold');
set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%1.1f'),'fontweight','bold');
ylim([0,200]);
h_lgd = legend(h,'15','16', '17', '18', '19', '20', '21','location','east');

% ==========================================================