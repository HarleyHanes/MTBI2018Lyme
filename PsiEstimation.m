function [eq,R,Psi] = PsiEstimation(params)
%PsiEstimation Calculates equilibrium values and R at various Psi(vaccination)
%              Values
%   Inputs: params- structure of parameters (psi and beta's calculated in
%   function)
%   Outputs: eq - structure with S,M,B elements for small, medium, big beta
%                values
%               - S,M,B each have a 10xc matrix where each row is an equilibrium
%                value for model compartment and each column is a psi
%                tested
%            R  - Vector of R values at each psi tested
%            Psi- Vector of psi values tested
BLarge=1.14;                            %Nearing 100% Nymphs infected w/out vaccination
BMed=.86;                               %W/out Vaccination ~65% Mice and ~85% Nymphs infected
BSmall=.68;                             %W/out Vaccination ~60% Mice and ~25% Nymphs Infected 
                                            %~Nymphal Prevalence seen in
                                            %Vaccine study
c=120;                                  %Number of psi values to test
Psi=linspace(0,12,c);                   %Create a vector of psi values to test
for k=1:c
    params.m.psi=Psi(k);                %Load psi value for interation
    for i=1:3                           %Run the model 3 times for each contact rate
        if i==1
            params.n.beta=BLarge;               %Load nymph contact rate
            params.l.beta=params.n.beta/.2;     %Adjust larvae contact rate
            params.m.beta=params.n.beta/.0968;  %Adjust mouse contact rate
            %Psi(k)=params.m.psi;
            sol=LymeSolutions(params);
            eq.B(:,k)=sol(:,end);
            R(i,k)=CalculateR(params);
            %R(i,k)=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(exp(-params.m.mu)+3*exp(-params.m.mu/4))/(exp(-3*params.m.mu/4)+3*exp(-params.m.mu));
            %R(i,k)=R(i,k)*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4);
        elseif i==2
            params.n.beta=BMed;
            params.l.beta=params.n.beta/.2;
            params.m.beta=params.n.beta/.0968;
            %Psi(k)=params.m.psi;
            sol=LymeSolutions(params);
            eq.M(:,k)=sol(:,end);
            R(i,k)=CalculateR(params);
            %R(i,k)=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu));
            %R(i,k)=R(i,k)*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4);
            %tau(k)=exp(-params.m.mu)+params.m.beta*params.l.beta/8*(1-exp(-params.m.mu))/(1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu));
        elseif i==3 
            params.n.beta=BSmall;
            params.l.beta=params.n.beta/.2;
            params.m.beta=params.n.beta/.0968;
            %Psi(k)=params.m.psi;
            sol=LymeSolutions(params);
            eq.S(:,k)=sol(:,end);
            R(i,k)=CalculateR(params);
            %R(i,k)=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu));
            %R(i,k)=R(i,k)*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4);
        end
    end
end
end

