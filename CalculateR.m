function R = CalculateR(params)
%CalculateR Calculates R (equilibrium staibility condition but not R0)
%Calculates using formulae calculated in paper
%   Inputs: params- structure of parameters
%   Output: R     -Equilibrium Stability Condtion
kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
b=exp(-params.m.mu);
R=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
end

