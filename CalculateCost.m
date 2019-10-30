function [Csaved] = CalculateCost(NI,Cbase,gamma,H_S)
%CalculateCost Calculates total cost of vaccination and Lyme cases 
%   Sums cost over the length of the solution and accounts for cost
%   increase due to inflation
%   Inputs: NI - vector of % infected nymphs at each time point
%                 ==(sol(10,:))
%           Cbase - cost without vaccination
%                 - input 0 to get a negative Cbase
%           gamma -vector of gammas to assess at
rho=.031;       % probability of human infection from contact
theta=3537.7;   % $ Cost of 1 case of Lyme Disease
%DALY=1.924;     % DALY cost of 1 case of Lyme Disease (DALY currently non-functional)
x=329.29;       % $ Cost per psi
R=1.017;        % Yearly cost increase

%Check input properties
if min(size(NI))~=1
    warning('NI is not a vector');
    keyboard
end
if size(Cbase)~=[1,1]
    warning('Cbase not a scaler');
end
if min(size(gamma))~=1
    warning('gamma is not a vector');
    keyboard
end

%Calculate Inflation
Inflation=R.^(tmax-1:-1:0); %we assume costs accrued at end of the year
Csaved=NaN(1,length(gamma));
for i=1:length(gamma)
    C=x*psi+theta*rho*gamma*H_S*NI(10,:);  %Calculate cost
    Cscaled=C.*Inflation;                  %Account for inflation
    Ctotal=sum(Cscaled);                   %Sum total cost over time
    Csaved(i)=Cbase-Ctotal;                  %Calculate change in cost
end

end


