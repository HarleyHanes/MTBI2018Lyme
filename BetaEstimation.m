function [eq,R0,Beta] = BetaEstimation(params,Beta_n_max,iter)
%BetaEstimation solves the system for variable values of Beta to estimate
%the accurate value
%Required Functions: LymeSolutions, R0finder
%Inputs
%   params- baselineparams inputed from baseline_params()
%   Beta_n_max- maximum value of Beta_n to be assessed
%   iter- number of different sets of Betas to solve infection for 
%Output
%   eq- 5 by iter matrix of equilibrium values for each compartment
%       (rows) at each assessed set of Beta values (columns)
%   R0-  iter length vector of R0 values for each set of Beta values,
%       calculated by solving our symbolic expresion for R0
%   Beta- 3 by iter matrix of Beta values assessed where rows are Beta_n,
%   Beta_l, then Beta_m


%    [c_l,c_n]=ContactRateCalc(params.m.popsize);
 %   c_m=c_n*(params.n.popsize/params.m.popsize); 
 
 
% Step 1: Comput Beta values to be assessed
     %linearly space Beta_n values from 0 to Beta_n_max
     Beta_n=linspace(0,Beta_n_max,iter);
     %define Beta_l to be 5 times Beta_n (Ratio of Beta_l to Beta_n 
     %found in literature values of Beta)
     Beta_l=Beta_n/.2;
     %define Beta_m to be 1/.0968 times Beta_n (Ratio of Beta_l to Beta_n 
     %found in literature values of Beta)
     Beta_m=Beta_n/.0968;
     %Organize into matrix
     Beta=[Beta_n; Beta_l; Beta_m];
%Step 2: Iterize
for k=1:iter
    %Step 2a: Substitute Parameters
    %Take the Beta values for the given iteration from the Beta matrix and
    %load them as the corresponding parameter values
    params.n.beta=Beta(1,k);
    params.l.beta=Beta(2,k);
    params.m.beta=Beta(3,k);
    %Step 2b: Solve System
    %Solve System using LymeSolutions
    sol=LymeSolutions(params);
    %Define equilibrium to be final LymeSolutions value
    eq(:,k)=sol(:,end);
    %Step 2c: Calculate R0
    R0(k)=R0finder(params);
end
end

