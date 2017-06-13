Lc = [0.015 0.0075 0.01 0.01 0.01];
d = [2.50E-03 2.50E-03 3.33E-03 (5/3000) 0.00025];
w = [0.025 0.025 0.025 0.025 0.025];
zeta = [5.274606666 4.922168903 4.834588912 5.62704443 51.47355896];
Lp = 0.005;

len=length(Lc);

x=cell(1,4);
for r=1:1:len
a = linspace(0,Lc(r),100);
x{r} = a;
end
beta_out = [0.1 0.1 0.1 0.1 0.1];
beta_in = [0.8 0.8 0.8 0.8 0.8];
rho = 1.225;
vo=1;

n = Lc./d;
Ac=d.^2.*pi./4;
Ai=d.*w;
Ae=d.*w;


% epsilon = zeros(len,1);
m=zeros(len,1);
% m_sq=zeros(len,1);
epsilon=(2-beta_in)./zeta./((Ai./Ac./n).^2);
m_sq=((2-beta_out)./(2-beta_in).*((Ai./Ae).^2)-1).*epsilon; 
m=sqrt(m_sq);
uc={};

for k=1:1:len
% pi_minus_pe{k}=((zeta/2).*((Ai./Ac./n(k)).^2)).*(((m(k)./sin(m(k))).^2).*(cos(m(k)*(1-x))-(epsilon_n(k)./(-m_sq(k))).*(cos(m(k)*x)-cos(m(k)*(1-x)))).^2);

pi_minus_pe{k}=((zeta(k)/2).*((Ai(k)/Ac(k)/n(k))^2)).*(((m(k)./sinh(m(k))).^2).*(cosh(m(k)*x{r}) + (epsilon(k)./m_sq(k)).*(cosh(m(k)*x{r})-cosh(m(k)*(1-x{r})))).^2).*(rho*vo^2);

% net_deltap_non_dim{k}=((zeta/2).*((Ai./Ac./n(k)).^2)).*(((m(k)./tanh(m(k))).^2) + epsilon(k).*(1 + epsilon(k)./(m_sq(k))).*((1-sech(m(k)))./tanh(m(k))).^2);
end
%% PLOTTING
plotStyle = {'b-','m-','g-'};
Lc_by_w = Lc./w;

% grid on;
% hold on;
% 
% figure;
% plot(x{1},pi_minus_pe{1},plotStyle{1},x{2},pi_minus_pe{2},plotStyle{2},...
%     x{3},pi_minus_pe{3},plotStyle{3},'LineWidth',1.2);



% legend(['L_c/w = ' num2str(Lc_by_w(1))],['L_c/w = ' num2str(Lc_by_w(2))], ...
%     ['L_c/w = ' num2str(Lc_by_w(3))],'Location','southwest','FontSize',9,'FontWeight','bold');
% xlabel('Linear location of channel \it (x)','FontSize',10,'FontWeight','bold')
% ylabel('\Delta P \it (x)','FontSize',10)

xcomp={}; pi_minus_pe_comp={};
xcomp{1} = [0.0012625 0.0037875 0.0063125 0.0088375 0.0113625 0.0139];
pi_minus_pe_comp{1} = [5.913092 6.092915 6.408052 6.86852 7.5590355 8.45639298];
%
xcomp{2} = [0 0.0013 0.0025 0.0038 0.0051 0.0063];
pi_minus_pe_comp{2} = [20.566652 21.057894 21.395907 22.207113 22.864066 23.9916748];
%
xcomp{3} = [0 0.00168165 0.0033633 0.00504495 0.006734925 0.0084249];
pi_minus_pe_comp{3} = [10.970925 11.367922 11.615892 12.266868 12.803279 14.6143353];
%
xcomp{4} = [0.00084335 0.00253005 0.00421675 0.00590345 0.00759015 0.0092852];
pi_minus_pe_comp{4} = [15.987436 16.232824 16.653054 17.254054 18.0385697 18.91306296];
%
xcomp{5} = [0.00441875 0.00467125 0.00492375 0.00921625 0.00946875 0.0097225];
pi_minus_pe_comp{5} = [257.9775 258.24779 259.25367 264.608561 265.390512 266.234428];

% figure(2)
% grid on;
% hold on;
% plot(x{2},pi_minus_pe{2},plotStyle{2});
% hold on;
% plot(xcomp{2},pi_minus_pe_comp{2},'r-o');
% % semilogy(xcomp{2},pi_minus_pe_comp{2},'b-o');    
%% YY TRY

figure(1);
hold on;
param1 = Lc./w;
param2 = d/Lp;

subplot(2,2,1)
hold on;
p=5;
 [AX,HA1,HA2] = plotyy(x{p},log(pi_minus_pe{p}),xcomp{p},log(pi_minus_pe_comp{p}),'plot');
% [AX,H1,H2] = plotyy(x{p},pi_minus_pe{p},xcomp{p},pi_minus_pe_comp{p},'plot');
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.2f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',9,'FontWeight','bold');
xlabel('$x/w$','Interpreter','Latex','FontSize',11,'FontWeight','bold');
grid on;
set(HA1,'LineStyle','--','LineWidth',2)
set(HA2,'LineStyle',':','marker','o','LineWidth',2)


subplot(222)
hold on;
p=3;
 [BX,HB1,HB2] = plotyy(x{p},log(pi_minus_pe{p}),xcomp{p},log(pi_minus_pe_comp{p}),'plot');
% [AX,H1,H2] = plotyy(x{p},pi_minus_pe{p},xcomp{p},pi_minus_pe_comp{p},'plot');
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.2f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',9,'FontWeight','bold');
xlabel('$x/w$','Interpreter','Latex','FontSize',11,'FontWeight','bold');
grid on;
set(HB1,'LineStyle','--','LineWidth',2)
set(HB2,'LineStyle',':','marker','o','LineWidth',2)


subplot(223)
hold on;
p=1;
 [CX,HC1,HC2] = plotyy(x{p},log(pi_minus_pe{p}),xcomp{p},log(pi_minus_pe_comp{p}),'plot');
% [AX,H1,H2] = plotyy(x{p},pi_minus_pe{p},xcomp{p},pi_minus_pe_comp{p},'plot');
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.2f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',9,'FontWeight','bold');
xlabel('$x/w$','Interpreter','Latex','FontSize',11,'FontWeight','bold');
grid on;
set(HC1,'LineStyle','--','LineWidth',2)
set(HC2,'LineStyle',':','marker','o','LineWidth',2)

%
hold on;
grid on;
%Ylabels
set(get(AX(1),'Ylabel'),'String','log(P_i - P_e) - Analytical','FontSize',8);
set(get(AX(2),'Ylabel'),'String','log(P_i - P_e) - Computational','FontSize',8);%'$Flow\ Unevenness={}^{{{u}_{MAX-Channel}}-{{u}_{HagenPoiseuille}}}/{}_{{{u}_{inlet}}}$', 'Interpreter', 'Latex')
%
% set(get(BX(1),'Ylabel'),'String','log(P_i - P_e) - Analytical','FontSize',9,'FontWeight','bold');
% set(get(BX(2),'Ylabel'),'String','log(P_i - P_e) - Computational','FontSize',9,'FontWeight','bold');%'$Flow\ Unevenness={}^{{{u}_{MAX-Channel}}-{{u}_{HagenPoiseuille}}}/{}_{{{u}_{inlet}}}$', 'Interpreter', 'Latex')
%
set(get(CX(1),'Ylabel'),'String','log(P_i - P_e) - Analytical','FontSize',8);
set(get(CX(2),'Ylabel'),'String','log(P_i - P_e) - Computational','FontSize',8);%'$Flow\ Unevenness={}^{{{u}_{MAX-Channel}}-{{u}_{HagenPoiseuille}}}/{}_{{{u}_{inlet}}}$', 'Interpreter', 'Latex')


%Use the line handles to set the LineStyle properties of the left- and right-side plots:
set(H1,'LineStyle','--','LineWidth',2)
set(H2,'LineStyle',':','marker','o','LineWidth',2)

%legend([H1;H2],'y1','y2');

% get current (active) axes property
scale = 0.041;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos);
