clc
close all
clear all

delta_T=43;
delta_P=50;

rho=1.225;
k=0.0242;
cp=1006.43;
mu=1.7894e-5;
alpha=k/rho/cp;

d = [0.0005 0.001 0.0025 0.005];
da= linspace(0.0005, 0.005, 1000);
Lp = [0.01 0.01 0.01 0.01];
for i=1:1:length(da)
    Lpa(i) = 0.01;
end

d_by_Lp = d./Lp;
da_by_Lpa = da./(Lpa);

Q1_pervol = rho*cp*delta_T*delta_P/32/mu.*da.^2./Lpa.^2;
R1 = delta_T./Q1_pervol;

Q2_pervol = 4*k*delta_T./sqrt(alpha*Lpa)./da.*(2*delta_P./rho)^0.25;
R2 = delta_T./Q2_pervol;

Qsim1_pervol = [1.10E+06 9.60E+05 1.76E+06 2.73E+05];
R_sim1 = delta_T./Qsim1_pervol;

Qsim2_pervol = [2.43E+06 1.99E+06 8.09E+05 4.17E+05];
R_sim2 = delta_T./Qsim2_pervol;

Qsim3_pervol = [5.17E+06 3.05E+06 1.65E+06 7.98E+05];
R_sim3 = delta_T./Qsim3_pervol;


grid on;
plot(da_by_Lpa,log(R1),da_by_Lpa,log(R2),d_by_Lp,log(R_sim1),d_by_Lp,log(R_sim2),d_by_Lp,log(R_sim3));
legend('d tends to 0','d tends to inf','Simulation Results');