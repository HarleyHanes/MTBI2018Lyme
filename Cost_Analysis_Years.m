function [C_vacc,C_case,C_total] = Cost_Analysis_Years(params,H_popsize,Beta,Psi,Time,exp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
params.m.psi=0;
params.n.beta=Beta(1)*365;
params.l.beta=Beta(2)*365;
params.m.beta=Beta(3)*365;
Init=LymeSolutions(params);
params.m.psi=Psi*365;
params.m.So=Init(1,end);
params.m.Io=Init(2,end);
params.n.So=Init(4,end);
params.n.Io=Init(5,end);
for i=1:length(Time)
    if Time(i)==0
        Infec(i)=Init(5,end);
    else
        params.tmax=Time(i);
        sol=LymeSolutions(params);
        Infec(i)=sol(5,end);
    end
    if exp==0
        C_vacc(i)=140556.4*Psi*Time(i);
    else
        C_vacc(i)=log(Psi)/log(.9999)*Time(i);
    end
end
gamma=.005*365/2;
rho=.031;
H=Infec.*gamma.*rho.*H_popsize./1000;
theta=5000;
C_case=theta*H;
C_total=C_case+C_vacc;
end

