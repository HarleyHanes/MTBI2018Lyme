function [C_vacc,C_case,DALY_time,Psi,I] = Cost_Analysis_Psi(params,H_popsize,gamma,year)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Conditions
params.tmax=year;
rho=.031;
theta=3537.7;   % $ Cost of 1 case of Lyme Disease (Needs to be rescaled for inflation)
DALY=1.924;     % DALY cost of 1 case of Lyme Disease (DALY currently non-functional)
x=329.29;

%Calculations
sol=LymeSolutions(params);
I_base=sol(5,end).*gamma.*rho.*H_popsize./1000;

[eq,R,Psi]=PsiEstimation(params);
infected=eq.M(5,:);
I=infected./1000*gamma.*rho.*H_popsize;
C_case=theta*I;
DALY_time=DALY*I;
C_vacc=x*Psi;

end

