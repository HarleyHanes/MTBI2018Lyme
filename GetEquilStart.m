function [params] = GetEquilStart(params,tol)
%GetEquilStart Gives parameter settings as equilibirum
%   Detailed explanation goes here
told=params.tmax;
params.tmax=1000;                
params.tol=tol;
[sol,~]=LymeSolutions(params);
if length(sol) < params.tmax
    warning("Lyme Solutions did not converge in $i iterations",params.tmax)
    keyboard
end
%Update parameters with 
params.m.So=sol(1,end);           %First row of sol plots susceptible mice
params.m.Io=sol(1,end);           %Second row of sol plots infected mice
params.m.Vo=sol(1,end);           %Third row of sol plots vaccinated mice
params.n.So=sol(1,end);           %Fourth row of sol plots susceptible nymphs
params.n.Io=sol(1,end);           %Fifth row of sol plots infected nymphs
params.tmax=told;               %Replace tmax with original

end

