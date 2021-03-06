% Note that Lp and w are held constant for these
Lc = [0.015 0.0075 0.01 20/3000 0.01];
d = [2.50E-03 2.50E-03 3.33E-03 (5/3000) 0.00025];
w = [0.025 0.025 0.025 0.025 0.025];
zeta = [5.274606666 4.922168903 4.834588912 5.128310069 51.47355896];
Lp = 0.005;
len=length(Lc);


param1 = Lc./w;
param2 = d/Lp;

x=cell(1,4);
for r=1:1:len
a = linspace(0,1,100);%(0,Lc(r),100);
x{r} = a;
end
beta_out = [0.1 0.1 0.1 0.1 0.1];
beta_in = [0.83 0.83 0.83 0.83 0.83];
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
for k=1:1:len
% pi_minus_pe{k}=((zeta/2).*((Ai./Ac./n(k)).^2)).*(((m(k)./sin(m(k))).^2).*(cos(m(k)*(1-x))-(epsilon_n(k)./(-m_sq(k))).*(cos(m(k)*x)-cos(m(k)*(1-x)))).^2);
vi{k} = 1 - (sinh(m(k)*x{r})./sinh(m(k))) + (epsilon(k)/m_sq(k))*(1-(sinh(m(k)*x{r})./sinh(m(k)))-(sinh(m(k)*(1-x{r}))./sinh(m(k))));
ve{k} = 1-vi{k};
end
%% PLOTTING
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
xcomp={}; vi_comp={}; xcomp_raw={}; ve_comp={};
%
xcomp_raw{1} = [-0.001275 0.0012625 0.0037875 0.0063125 0.0088375 0.0113625 0.0139];
vi_comp{1} = [0.984634500000000,0.881735000000000,0.754877800000000,0.593874500000000,0.404017500000000,0.212452500000000,0.0199841700000000];
ve_comp{1} = [0.000339652000000000,0.0885920700000000,0.224447000000000,0.384436900000000,0.569579500000000,0.767458500000000,0.997379400000000];
%
xcomp_raw{2} = [-0.001275 0 0.0013 0.0025 0.0038 0.0051 0.0063];
vi_comp{2} = [0.986207000000000,0.857502200000000,0.722725800000000,0.570990500000000,0.376892200000000,0.191230800000000,0.0508632500000000];
ve_comp{2} = [0.000333454000000000,0.103322200000000,0.246303400000000,0.417214000000000,0.589619100000000,0.786550200000000,0.993029700000000];
%
xcomp_raw{3} = [-0.0016983 0 0.00168165 0.0033633 0.00504495 0.006734925 0.0084249];
vi_comp{3} = [0.996366500000000,0.871731900000000,0.735737800000000,0.579420000000000,0.395846000000000,0.201104700000000,0.0626588500000000];
ve_comp{3} = [-0.000408838000000000,0.0898148800000000,0.224893200000000,0.396329500000000,0.566168500000000,0.764357300000000,1.08652500000000];
%
xcomp_raw{4} = [-0.0008517 0.00084335 0.00253005 0.00421675 0.00590345 0.00759015 0.0092852];
vi_comp{4} = [0.979631800000000,0.870518100000000,0.725240900000000,0.577383700000000,0.388141000000000,0.192991200000000,0.0140873400000000];
ve_comp{4} = [0.00110600600000000,0.109404500000000,0.253950400000000,0.416624300000000,0.603208700000000,0.796686500000000,0.996123900000000];
%
xcomp_raw{5} = [-0.000127500000000000,0.000126250000000000,0.000378750000000000,0.000631250000000000,0.000883750000000000,0.00441875000000000,0.00467125000000000,0.00492375000000000,0.00921625000000000,0.00946875000000000,0.00972250000000000];
vi_comp{5} = [0.999990000000000,0.927161700000000,0.930187800000000,0.903716300000000,0.874733000000000,0.559445100000000,0.535288200000000,0.506682500000000,0.0912705800000000,0.0632149300000000,0.0365957000000000];
ve_comp{5} = [-0.000421215000000000,0.0212934300000000,0.0443775300000000,0.0678768700000000,0.0895919300000000,0.415045000000000,0.437709300000000,0.460188600000000,0.882495200000000,0.904280800000000,0.922254500000000];
for t=1:1:len
xcomp{t} = (xcomp_raw{t}-min(xcomp_raw{t}))./(max(xcomp_raw{t})-min(xcomp_raw{t}));
end
%% SINGLE PLOT
p=2;
figure(p)
grid on;
hold on;
%vi
plot(x{p},vi{p},'b--','LineWidth',1.2);
hold on;
plot(xcomp{p},vi_comp{p},'r:','marker','o','LineWidth',1.7);
%ve
plot(x{p},ve{p},'k--','LineWidth',1.2);
hold on;
plot(xcomp{p},ve_comp{p},'dg:','marker','o','LineWidth',1.7);
% Labelling...
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.3f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
xlabel('$x^*$','Interpreter','Latex','FontSize',12,'FontWeight','bold');
legend('v_i Analytic','v_i Comp.','v_e Analytic','v_e Comp.')
p=2;
figure(p)
grid on;
hold on;
%vi
plot(x{p},vi{p},'b--','LineWidth',1.2);
hold on;
plot(xcomp{p},vi_comp{p},'r:','marker','o','LineWidth',1.7);
%ve
plot(x{p},ve{p},'k--','LineWidth',1.2);
hold on;
plot(xcomp{p},ve_comp{p},'dg:','marker','o','LineWidth',1.7);
% Labelling...
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.3f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
xlabel('$x^*$','Interpreter','Latex','FontSize',14,'FontWeight','bold');
ylabel('Non-dim. velocity','Interpreter','Latex','FontSize',14,'FontWeight','bold');
legend('v_i Analytic','v_i Comp.','v_e Analytic','v_e Comp.')
%% SINGLE PLOT
p=4;
figure(p)
grid on;
hold on;
%vi
plot(x{p},vi{p},'b--','LineWidth',1.2);
hold on;
plot(xcomp{p},vi_comp{p},'r:','marker','o','LineWidth',1.7);
%ve
plot(x{p},ve{p},'k--','LineWidth',1.2);
hold on;
plot(xcomp{p},ve_comp{p},'dg:','marker','o','LineWidth',1.7);
% Labelling...
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.3f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
xlabel('$x^*$','Interpreter','Latex','FontSize',14,'FontWeight','bold');
ylabel('Non-dim. velocity','Interpreter','Latex','FontSize',14,'FontWeight','bold');
legend('v_i Analytic','v_i Comp.','v_e Analytic','v_e Comp.')
%% SINGLE PLOT
p=5;
figure(p)
grid on;
hold on;
%vi
plot(x{p},vi{p},'b--','LineWidth',1.2);
hold on;
plot(xcomp{p},vi_comp{p},'r:','marker','o','LineWidth',1.7);
%ve
plot(x{p},ve{p},'k--','LineWidth',1.2);
hold on;
plot(xcomp{p},ve_comp{p},'dg:','marker','o','LineWidth',1.7);
% Labelling...
obj=title(sprintf('$m^2 = %0.4f$; $L_c/w = %0.3f$; $d/L_p = %0.3f$',m_sq(p),param1(p),param2(p)));
set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
xlabel('$x^*$','Interpreter','Latex','FontSize',14,'FontWeight','bold');
ylabel('Non-dim. velocity','Interpreter','Latex','FontSize',14,'FontWeight','bold');
legend('v_i Analytic','v_i Comp.','v_e Analytic','v_e Comp.')