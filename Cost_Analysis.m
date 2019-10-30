function Csaved = Cost_Analysis(params,H_popsize,gamma,year,psi)
%Cost_Analysis Estimates the cost of vaccinations and lyme disease
%treatment
%   Inputs: params - structure of parameters
%           H_popsize - Human risk size
%           year - After how many years do we assess cost?
%           exp- bolean for whether to calculate C_vac with exponential
%   Outputs: C_vacc - cost of vaccination
%            C_case - cost of Lyme disease
%            Psi-

%Get base conditions
params.psi=0;                           %Set vaccination equal to 0
%Get equilibrium starting populations
    params=GetEquilStart(params,10^(-8));   %Set params to equil condition
%Calculate cost without vaccination
    params.tmax=year;
    [solbase,~]=LymeSolutions(params,0);            %Calculate solution without stopping at convervence
    if norm(solbase(6:10,1)-solbase(6:10,end)) > 10^(-7) %Check Equilibrium
        warning("Equil parameters not at equilibrium")
        keyboard
    end
    if length(solbase(1,:))~=year                        %Check correct # years
        warning("Base case solution length %i, expected %i", length(solbase(1,:)),year)
        keyboard
    end
    Cbase = -CalculateCost(solbase(10,:),0,gamma,H_popsize);       %Calculate base cost
    if Cbase < 0
        warning("Base cost is negative")
    end
    
    
%Initialize Csaved
Csaved=NaN(length(psi),length(gamma)); 
for ipsi=1:length(psi)
    params.psi=psi(ipsi);
    [sol,~]=LymeSolution(params,0);
    Csaved(ipsi,:)=CalculateCost(sol(10,:),Cbase,gamma,H_popsize);
end
end

