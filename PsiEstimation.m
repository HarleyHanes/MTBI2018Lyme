function [eq,R,Psi] = PsiEstimation(params)
%PsiEstimation Summary of this function goes here
%   Detailed explanation goes here
BLarge=1.14;                            %Nearing 100% Nymphs infected w/out vaccination
BMed=.86;                               %W/out Vaccination ~65% Mice and ~85% Nymphs infected
BSmall=.68;                             %W/out Vaccination ~60% Mice and ~25% Nymphs Infected 
                                        %~Nymphal Prevalence seen in
                                        %Vaccine study
c=200;
P=linspace(0,12,c);
%R=.5*((1-exp(-params.m.mu))*(exp(-params.m.mu)+3*exp(-params.m.mu/4))*params.l.beta*params.m.beta*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4)/(8*(3*exp(-params.m.mu)+exp(-params.m.mu+exp(-3*params.m.mu/4)))*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4))))
for k=1:c
    %j=(k-1)*.4;
    %params.m.psi=.002283416825*365*j;
    params.m.psi=P(k);
    for i=1:3
        if i==1
            params.n.beta=BLarge;
            params.l.beta=params.n.beta/.2;
            params.m.beta=params.n.beta/.0968;
            Psi(k)=params.m.psi;
            sol=LymeSolutions(params);
            eq.B(:,k)=sol(:,end);
            kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
            a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
            b=exp(-params.m.mu);
            R(i,k)=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            %R(i,k)=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(exp(-params.m.mu)+3*exp(-params.m.mu/4))/(exp(-3*params.m.mu/4)+3*exp(-params.m.mu));
            %R(i,k)=R(i,k)*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4);
        elseif i==2
            params.n.beta=BMed;
            params.l.beta=params.n.beta/.2;
            params.m.beta=params.n.beta/.0968;
            Psi(k)=params.m.psi;
            sol=LymeSolutions(params);
            eq.M(:,k)=sol(:,end);
            kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
            a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
            b=exp(-params.m.mu);
            R(i,k)=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            %R(i,k)=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu));
            %R(i,k)=R(i,k)*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4);
            %tau(k)=exp(-params.m.mu)+params.m.beta*params.l.beta/8*(1-exp(-params.m.mu))/(1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu));
        elseif i==3 
            params.n.beta=BSmall;
            params.l.beta=params.n.beta/.2;
            params.m.beta=params.n.beta/.0968;
            Psi(k)=params.m.psi;
            sol=LymeSolutions(params);
            eq.S(:,k)=sol(:,end);
            kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
            a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
            b=exp(-params.m.mu);
            R(i,k)=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            %R(i,k)=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu));
            %R(i,k)=R(i,k)*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4);
        end
    end
end

end

