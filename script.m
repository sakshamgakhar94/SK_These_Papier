kmax=5;
for k=1:kmax
    %% Specification of the porous fin (all lenths in SI)
    layers=1;%%------------DONOT MAKE A DIFFERENCE--------------%%.
    s=linspace(0,1,5);
    
    %k=1;%%-------TBV-------------TBV-------------TBV-----------------%%.
    beta_i=s(k) %%-------TBV-------------TBV-------------TBV-----------------%
    beta_e=s(kmax-k+1) %%-------TBV-------------TBV-------------TBV---------%%.
    
    zeta=400;%%-------TBV-------------TBV-------------TBV---------%%.
    
    ppi=20;         %%-------TBV-------------TBV-------------TBV---------%%.
    L=0.04;%width of the channel into plane (fixed)
    n=L*ppi/0.0254;% note that L/n=0.0254/ppi
    n_actual=floor(n);% no of cylinders in 1 row
    x=linspace(0,1,50);
    
    %% Imposing constraint on conduit width a & Dia of  cyl pores: ratio: alpha
    alpha=40;   %%-------TBV-------------TBV-------------TBV---------%%.
    a=0.0254/ppi*alpha;
    %note that height of the channel = layers*L/n
    
    %% Areas of the channel (Ac) and the inlet and outelt conduits (Ai, Ae)
    Ac=((L/n)^2)*pi/4;
    Ai=(a/2)*((L/n)*layers);%b=(L/n)*layers
    Ae=(a/2)*((L/n)*layers);%b=(L/n)*layers %% FOR NOW SAME AS Ai
    
    %% Define parameters epsilon and m_sq
    epsilon=(2-beta_i)/zeta/((Ai/Ac/n)^2);
    m_sq=((2-beta_e)/(2-beta_i)*((Ai/Ae)^2)-1)*epsilon;
    
    %% Formulation for 2 cases of m>0 and m<0 (note quantities are all non dim)
    if m_sq>0
        m=sqrt(m_sq);
        %%Axial velocities of the conduit at intake and outlet of the channels
        wi=1-(sinh(m*x)/sinh(m))+(epsilon/m_sq)*(1-(sinh(m*x)/sinh(m))-(sinh(m*(1-x))/sinh(m)));
        we=Ai/Ae*(1-wi);
        %%Channel flow distribution
        uc=(Ai/Ac/n)*((m/sinh(m))*(cosh(m*x)+(epsilon/m_sq)*(cosh(m*x)-cosh(m*(1-x)))));
        % pressure drop in the channels
        pi_minus_pe=((zeta/2)*((Ai/Ac/n)^2))*(((m/sinh(m))^2)*(cosh(m*x)+(epsilon/m_sq)*(cosh(m*x)-cosh(m*(1-x)))).^2);
        %%Total pressure drop
        net_deltap_non_dim=((zeta/2)*((Ai/Ac/n)^2))*(((m/tanh(m))^2)+epsilon*(1+epsilon/(m_sq))*((1-sech(m))/tanh(m))^2);
    elseif m_sq<0
        m=sqrt(-m_sq);% we dnot introduce the m_prime
        epsilon_n=(1-beta_e)*((Ai/Ae)^2)/(zeta*(Ai/Ac/n)^2);%this is the epsilon prime
        %%Axial velocities of the conduit at intake and outlet of the channels
        wi=(sin(m*(1-x))/sin(m))-(epsilon_n/(-m_sq))*(1-(sin(m*x)/sin(m))-(sin(m*(1-x))/sin(m)));
        we=Ai/Ae*(1-wi);
        %%Channel flow distribution
        uc=(Ai/Ac/n)*((m/sin(m))*(cos(m*(1-x))-(epsilon_n/(-m_sq))*(cos(m*x)-cos(m*(1-x)))));
        % pressure drop in the channels
        pi_minus_pe=((zeta/2)*((Ai/Ac/n)^2))*(((m/sin(m))^2)*(cos(m*(1-x))-(epsilon_n/(-m_sq))*(cos(m*x)-cos(m*(1-x)))).^2);
        %%Total pressure drop
        net_deltap_non_dim=((zeta/2)*((Ai/Ac/n)^2))*(((m/tan(m))^2)+epsilon*epsilon_n/(-m_sq)*((1-sec(m))/tan(m))^2);
    elseif m==0
        %%Axial velocities of the conduit
        wi=1-x+epsilon/2*x*(1-x);
        we=(Ai/Ae)*(1-wi);
        %%Channel flow distribution
        uc=(Ai/Ac/n)*(1-epsilon/2*(1-2*x));
        % pressure drop in the channels
        pi_minus_pe=((zeta/2)*((Ai/Ac/n)^2))*((1-epsilon/2*(1-2*x))^2);
        %%Total pressure drop
        net_deltap_non_dim=((zeta/2)*((Ai/Ac/n)^2))*(1+(epsilon^2)/4);
    end
    
    %NET_PRESS_DROP=[];%generate the table at a fixed ppi, layers, alpha, zeta
    NET_PRESS_DROP(k,1)=net_deltap_non_dim;
    NET_PRESS_DROP(k,2)=beta_i;% beta_i
    NET_PRESS_DROP(k,3)=beta_e;% beta_e
    NET_PRESS_DROP(k,4)=zeta;% beta_e
    NET_PRESS_DROP(k,5)=layers;% beta_e
    NET_PRESS_DROP(k,6)=alpha;% beta_e
    NET_PRESS_DROP(k,7)=ppi;% beta_e
    NET_PRESS_DROP(k,8)=m_sq;% beta_e
    NET_PRESS_DROP(k,9)=max(uc);% beta_e
    NET_PRESS_DROP(k,10)=min(uc);% beta_e
    NET_PRESS_DROP(k,11)=mean(uc);% beta_e
    NET_PRESS_DROP(k,12)=(max(uc)-min(uc))*100/mean(uc);% beta_e
    
%     %% Plotting ALL PLOTS ARE NON-DIMENSIONAL
%     plotStyle = {'b-','k:','r.','g*','c^'}; % add as many as you need
%     hold on;
%     hold on;
%     figure(1);
%     plot(x,uc,plotStyle{k});
%     hold on;
%     xlabel({'X (non-dimensional)'},'FontSize',10,'FontWeight','bold');
%     ylabel({'uc (Channel velocity)'},'FontSize',10,'FontWeight','bold');
%     legendinfo{k}=[sprintf('\\beta _i = %0.2f; \\beta _e = %0.2f \n ',beta_i,beta_e)];%'Location','northwest','FontSize',10,'FontWeight','bold'
%     legend(legendinfo,'Location','northwest','FontSize',10,'FontWeight','bold')
%     obj=title(sprintf('$ppi= %d$; $layers= %d$; $\\displaystyle \\frac{Channel Width}{Pore Diameter}= %d$; $\\zeta= %d$',ppi,layers,alpha,zeta));%;
%     set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
%     hold off;
%     
%     %%--------------------INLET Conduit velocity-----------------------------------%%
%     
%     hold on;
%     figure(2);
%     hold on;
%     plot(x,wi,plotStyle{k});
%     hold on;
%     hold on;
%     xlabel({'X (non-dimensional)'},'FontSize',10,'FontWeight','bold');
%     ylabel({'w_{i}(Inlet Conduit velocity)'},'FontSize',10,'FontWeight','bold');
%     legendinfo{k}=[sprintf('\\beta _i = %0.2f; \\beta _e = %0.2f \n ',beta_i,beta_e)];%'Location','northwest','FontSize',10,'FontWeight','bold'
%     legend(legendinfo,'Location','northeast','FontSize',10,'FontWeight','bold')
%     obj=title(sprintf('$ppi= %d$; $layers= %d$; $\\displaystyle \\frac{Channel Width}{Pore Diameter}= %d$; $\\zeta= %d$',ppi,layers,alpha,zeta));%;
%     set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
%     hold off;
%     
%     %%--------------------EXIT OUTL:UET Conduit velocity-----------------------------------%%
%     
%     hold on;
%     figure(3);
%     plot(x,we,plotStyle{k});
%     hold on;
%     hold on;
%     xlabel({'X (non-dimensional)'},'FontSize',10,'FontWeight','bold');
%     ylabel({'w_{e}(Exit Conduit velocity)'},'FontSize',10,'FontWeight','bold');
%     legendinfo{k}=[sprintf('\\beta _i = %0.2f; \\beta _e = %0.2f \n ',beta_i,beta_e)];%'Location','northwest','FontSize',10,'FontWeight','bold'
%     legend(legendinfo,'Location','northwest','FontSize',10,'FontWeight','bold')
%     obj=title(sprintf('$ppi= %d$; $layers= %d$; $\\displaystyle \\frac{Channel Width}{Pore Diameter}= %d$; $\\zeta= %d$',ppi,layers,alpha,zeta));%;
%     set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
%     hold off;
%     
%     %%--------------------loacal pressure drops-----------------------------------%%
%     
%     hold on;
%     figure(4);
%     plot(x,pi_minus_pe,plotStyle{k});
%     hold on;
%     hold on;
%     xlabel({'X (non-dimensional)'},'FontSize',10,'FontWeight','bold');
%     ylabel({'P_{i} - P_{e} (Local pressure drop)'},'FontSize',10,'FontWeight','bold');
%     legendinfo{k}=[sprintf('\\beta _i = %0.2f; \\beta _e = %0.2f \n ',beta_i,beta_e)];%'Location','northwest','FontSize',10,'FontWeight','bold'
%     legend(legendinfo,'Location','northwest','FontSize',10,'FontWeight','bold')
%     obj=title(sprintf('$ppi= %d$; $layers= %d$; $\\displaystyle \\frac{Channel Width}{Pore Diameter}= %d$; $\\zeta= %d$',ppi,layers,alpha,zeta));%;
%     set(obj,'Interpreter','Latex','FontSize',12,'FontWeight','bold');
%     hold off;
end


%% Displaying data as a table: not used
plot_for_beta_var=0;
if plot_for_beta_var==1
    colheadings = {'Sr. No.','NET_PRESSURE_DROP','\beta_{i}','\beta_{e}'};
    rowheadings = {'','','', '',''} ;
    data=zeros(kmax,4);
    for i=1:1:kmax
        data(i,1)=i ;
        data(i,2)=NET_PRESS_DROP(i,1);
        data(i,3)=NET_PRESS_DROP(i,2);;
        data(i,4)=NET_PRESS_DROP(i,3);;
    end
    
    %To format the first number in each row as a decimal (%d), the second number %16.4f, and the third as %16.5E do the following:
    
    wid = 17;
    fms = {'i','.8f','.4f','.4f'};
    %In this case 16 will be the field width, and '.5E' is what to use for the fms argument
    fileID = 1;
    displaytable(data,colheadings,wid,fms,rowheadings,fileID);
end