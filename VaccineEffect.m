function [eq,R,Psi] = VaccineEffect(params,years)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
BLarge=1.14;                              %Nearing 100% Nymphs infected w/out vaccination
BMed=.86;                                %W/out Vaccination ~65% Mice and ~85% Nymphs infected
BSmall=.68;                               %W/out Vaccination ~60% Mice and ~25% Nymphs Infected 
                                        %~Nymphal Prevalence seen in
                                        %Vaccine study
c=100;
%calculating endemic equilibiri and initial conditions
params.m.psi=0;
params.n.beta=BSmall;
params.l.beta=params.n.beta/.2;
params.m.beta=params.n.beta/.0968;
sol=LymeSolutions(params);
eq.S.nymph(1,1:c)=sol(10,end);
eq.S.mouse(1,1:c)=sol(7,end);
eq.S.vacc(1,1:c)=sol(8,end);
MIo.small=sol(2,end);
NIo.small=sol(5,end);

params.n.beta=BMed;
params.l.beta=params.n.beta/.2;
params.m.beta=params.n.beta/.0968;
sol=LymeSolutions(params);
eq.M.nymph(1,1:c)=sol(10,end);
eq.M.mouse(1,1:c)=sol(7,end);
eq.M.vacc(1,1:c)=sol(8,end);
MIo.med=sol(2,end);
NIo.med=sol(5,end);

params.n.beta=BLarge;
params.l.beta=params.n.beta/.2;
params.m.beta=params.n.beta/.0968;
sol=LymeSolutions(params);
eq.B.nymph(1,1:c)=sol(10,end);
eq.B.mouse(1,1:c)=sol(7,end);
eq.B.vacc(1,1:c)=sol(8,end);
MIo.big=sol(2,end);
NIo.big=sol(5,end);

% Determining Psi

Psi=linspace(0,15,c);

for k=1:c
    params.m.psi=Psi(k);
    for j=[2 3 4]
        if j==2
            params.tmax=years(1);
        elseif j==3
            params.tmax=years(2);
        elseif j==4
            params.tmax=years(3);
        end
        for i=1:3
            if i==1
                params.n.beta=BLarge;
                params.l.beta=params.n.beta/.2;
                params.m.beta=params.n.beta/.0968;
                params.m.Io=MIo.big;
                params.n.Io=NIo.big;
                sol=LymeSolutions(params);
                eq.B.nymph(j,k)=sol(10,end);
                eq.B.mouse(j,k)=sol(7,end);
                eq.B.vacc(j,k)=sol(8,end);
                kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
                a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
                b=exp(-params.m.mu);
                R(i,k)=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            elseif i==2
                params.n.beta=BMed;
                params.l.beta=params.n.beta/.2;
                params.m.beta=params.n.beta/.0968;
                params.m.Io=MIo.med;
                params.n.Io=NIo.med;
                sol=LymeSolutions(params);
                eq.M.nymph(j,k)=sol(10,end);
                eq.M.mouse(j,k)=sol(7,end);
                eq.M.vacc(j,k)=sol(8,end);
                kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
                a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
                b=exp(-params.m.mu);
                R(i,k)=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            elseif i==3 
                params.n.beta=BSmall;
                params.l.beta=params.n.beta/.2;
                params.m.beta=params.n.beta/.0968;
                params.m.Io=MIo.small;
                params.n.Io=NIo.small;
                sol=LymeSolutions(params);
                eq.S.nymph(j,k)=sol(10,end);
                eq.S.mouse(j,k)=sol(7,end);
                eq.S.vacc(j,k)=sol(8,end);
                kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
                a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
                b=exp(-params.m.mu);
                R(i,k)=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            end
        end
    end
end

end

