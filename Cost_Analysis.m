function [C_vacc,C_case,Psi] = Cost_Analysis(params,H_popsize,year,exp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
params.tmax=year;
[eq,R,Psi]=PsiEstimation(params);
infected=eq.S(5,:);
gamma=.005*365;
rho=.031;
H=infected.*gamma.*rho.*H_popsize./1000;
theta=5000;
C_case=theta*H;
if exp==0
    C_vacc=140556.4*Psi*year;
else
    C_vacc=log(Psi)/log(.9999)*year;
end
C_total=C_case+C_vacc;
end

