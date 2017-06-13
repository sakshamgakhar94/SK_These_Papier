clc
close all
% clear all

%deltas
delta_T=43;
delta_P=50;

%FLUID PROPS
rho=1.225;
k=0.0242;
cp=1006.43;
mu=1.7894e-5;
alpha=k/rho/cp;

%GEOMETRY PARAMS
w = [0.025 0.025 0.025];
Lc = [0.025 0.01 0.004];
d = [0.0005 0.001 0.0025 0.005];
da= linspace(0.0002, 0.005, 1000);
Lp = [0.01 0.01 0.01 0.01];
for i=1:1:length(da)
    Lpa(i) = 0.01;
end

dbyLp = d./Lp;
da_by_Lpa = da./(Lpa);

d_by_Lp={};
d_by_Lp{1} = dbyLp;
d_by_Lp{2} = [0.0002 0.0005 0.001 0.0025 0.005]./[0.01 0.01 0.01 0.01 0.01];
d_by_Lp{3} = dbyLp;
%ANALYTICAL Q AND R : 2 LIMITS

% 0 Limit
Q1_pervol = rho*cp*delta_T*delta_P/32/mu.*da.^2./Lpa.^2;
R1 = delta_T./Q1_pervol;

% INF Limit
Q2_pervol = 4*k*delta_T./sqrt(alpha*Lpa)./da.*(2*delta_P./rho)^0.25;
R2 = delta_T./Q2_pervol;

% Simulation Results
Qsim_pervol={};
unevenness = {};
R_sim={};

Qsim_pervol{1} = [1.10E+06 14.330E+05 1.76E+06 2.73E+05]; %% 2nd entry was 9.6e5
R_sim{1} = delta_T./Qsim_pervol{1};
% unevenness{1} = [0.404420911 0.236120911 0 -0.061782724];

Qsim_pervol{2} = [2.50E+06 2.43E+06 1.99E+06 8.09E+05 4.17E+05];
R_sim{2} = delta_T./Qsim_pervol{2};
% unevenness{2} = [-0.65 -0.784097724 -1.008397724 -1.166968096 -1.713897724];%% 1st entry is arbit

Qsim_pervol{3} = [5.17E+06 3.05E+06 1.65E+06 7.98E+05];
R_sim{3} = delta_T./Qsim_pervol{3};
% unevenness{3} = [-1.260492487 -2.09053963 -4.125494309 -5.004794309];
%% PLOTTING
plotStyle = {'k--','r--','b-o','m-s','g-^'};
grid on;
hold on;

plot(da_by_Lpa,log(R1),plotStyle{1},da_by_Lpa,log(R2),plotStyle{2},'LineWidth',2);

plot(d_by_Lp{1},log(R_sim{1}),plotStyle{3},d_by_Lp{2},log(R_sim{2}),plotStyle{4},...
    d_by_Lp{3},log(R_sim{3}),plotStyle{5},'LineWidth',1.2);

Lc_by_w = Lc./w;

legend('d/L_p\rightarrow 0 limit','d/L_p\rightarrow \infty limit',['L_c/w = ' num2str(Lc_by_w(1))],...
    ['L_c/w = ' num2str(Lc_by_w(2))],['L_c/w = ' num2str(Lc_by_w(3))],'Location','southwest','FontSize',9,'FontWeight','bold');
xlabel('d/L_p','FontSize',10,'FontWeight','bold')
ylabel('log(R_{thermal})  [K-m^{3}/W]','FontSize',10)
 
%% YY PLOTTING
% p = 3 ; % Case # = p
% figure(1);
% [AX,H1,H2] = plotyy(d_by_Lp{p},R_sim{p},d_by_Lp{p},unevenness{p},'plot');
% 
% %You can use the handles returned by plotyy to label the axes and set the line styles used for plotting.
% % With the axes handles you can specify the YLabel properties of the left- and right-side y-axis:
% hold on;
% grid on;
% %Ylabels
% set(get(AX(1),'Ylabel'),'String','R_{thermal} (K-m^{3}/W)');
% set(get(AX(2),'Ylabel'),'String','Flow Unevenness');%'$Flow\ Unevenness={}^{{{u}_{MAX-Channel}}-{{u}_{HagenPoiseuille}}}/{}_{{{u}_{inlet}}}$', 'Interpreter', 'Latex')
% 
% %Xlabel
% % xlab='d/L_{p}';
% % obj_label=xlabel(xlab);
% % set(obj_label,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
% xlabel('$d/L_{p}$','Interpreter','Latex','FontSize',12,'FontWeight','bold');

% %TiTle
% param = Lc./w;
% obj=title(sprintf('$Lc/w = %0.3f$;',param(p)));
% set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
% 
% %Use the line handles to set the LineStyle properties of the left- and right-side plots:
% set(H1,'LineStyle',':','marker','s','LineWidth',2)
% set(H2,'LineStyle',':','marker','o','LineWidth',2)
% 
% %legend([H1;H2],'y1','y2');
% 
% % get current (active) axes property
% scale = 0.041;
% pos = get(gca, 'Position');
% pos(2) = pos(2)+scale*pos(4);
% pos(4) = (1-scale)*pos(4);
% set(gca, 'Position', pos)
% 
% %
% % fname = ['R_uneven_' num2str(p) '.fig'];
% % savefig(fname)
% % ename = ['R_uneven_' num2str(p) '.eps'];
% % saveas(gcf,ename)
% % jname = ['R_uneven_' num2str(p) '.jpg'];
% % saveas(gcf,jname)
