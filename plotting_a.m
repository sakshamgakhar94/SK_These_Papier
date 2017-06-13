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
w = [0.005 0.015 0.05];
Lc = [0.015 0.015 0.015];
d = [0.0005 0.001 0.0025 0.005];
da= linspace(0.0003, 0.005, 1000);
Lp = [0.01 0.01 0.01 0.01];
for i=1:1:length(da)
    Lpa(i) = 0.01;
end

dbyLp = d./Lp;
da_by_Lpa = da./(Lpa);

d_by_Lp={};
d_by_Lp{1} = [0.0003 0.0005 0.001 0.0025 0.005]./[0.01 0.01 0.01 0.01 0.01];
d_by_Lp{2} = dbyLp;
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

Qsim_pervol{1} = [1.94E+06 1.94E+06 1.49E+06 9.58E+05 5.64E+05]; %% 2nd entry was 9.6e5
R_sim{1} = delta_T./Qsim_pervol{1};
unevenness{1} = [0.085 0.161345868 0.359109154 0.565205658 0.372173637];

Qsim_pervol{2} = [1.73E+06 1.49E+06 7.82E+05 4.15E+05];
R_sim{2} = delta_T./Qsim_pervol{2};
unevenness{2} = [0.11642391 0.183890682 0.158531625 -0.202479089];

Qsim_pervol{3} = [1.63E+06 1.11E+06 4.99E+05 2.52E+05];
R_sim{3} = delta_T./Qsim_pervol{3};
unevenness{3} = [-0.195347324 -0.451598127 -1.4 -2.378263632];

% grid on;
% plot(da_by_Lpa,log(R1),da_by_Lpa,log(R2),d_by_Lp,log(R_sim{1}),d_by_Lp,log(R_sim{2}),d_by_Lp,log(R_sim{3}));
% legend('d tends to 0','d tends to inf','Simulation Results');


%% PLOTTING
p = 3 ; % Case # = p
figure(p);
[AX,H1,H2] = plotyy(d_by_Lp{p},R_sim{p},d_by_Lp{p},unevenness{p},'plot');

%You can use the handles returned by plotyy to label the axes and set the line styles used for plotting.
% With the axes handles you can specify the YLabel properties of the left- and right-side y-axis:
hold on;
grid on;
%Ylabels
set(get(AX(1),'Ylabel'),'String','R_{thermal} [K-m^{3}/W]');
set(get(AX(2),'Ylabel'),'String','Flow Unevenness');%'$Flow\ Unevenness={}^{{{u}_{MAX-Channel}}-{{u}_{HagenPoiseuille}}}/{}_{{{u}_{inlet}}}$', 'Interpreter', 'Latex')

%Xlabel
% xlab='d/L_{p}';
% obj_label=xlabel(xlab);
% set(obj_label,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
xlabel('$d/L_{p}$','Interpreter','Latex','FontSize',12,'FontWeight','bold');

%TiTle
param = w./Lc;
obj=title(sprintf('$w/L_c = %0.3f$;',param(p)));
set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');

%Use the line handles to set the LineStyle properties of the left- and right-side plots:
set(H1,'LineStyle',':','marker','s','LineWidth',2)
set(H2,'LineStyle',':','marker','o','LineWidth',2)

%legend([H1;H2],'y1','y2');

% get current (active) axes property
scale = 0.041;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

%
% fname = ['R_uneven_' num2str(p) '.fig'];
% savefig(fname)
% ename = ['R_uneven_' num2str(p) '.eps'];
% saveas(gcf,ename)
% jname = ['R_uneven_' num2str(p) '.jpg'];
% saveas(gcf,jname)
