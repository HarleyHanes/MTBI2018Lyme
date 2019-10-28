function [params] = Lymeparams()
%Lymeparams Enters all parameter values for numerical simulations
%   No computation is done in this function, only holds parameter, initial
%   conditions, and time values
%% Intial conditions
params.l.popsize=0;
params.n.popsize=1000;
params.a.popsize=0;
params.m.popsize=50;
params.m.Io=10;
params.m.So=params.m.popsize-params.m.Io;
params.m.Vo=0;
params.n.Io=100;
params.n.So=params.n.popsize-params.n.Io;
%% Mouse Parameters
% All rate parameters are in units 1/year
params.m.mu=.012*365;                                %Source:
params.m.lambda=4*params.m.popsize*(1-exp(-params.m.mu))/(exp(-3*params.m.mu/4)+3);                             %Source:
params.m.psi=0;%.002283416825*365;                         %Source: Field Trial, Captures/ Nights of trap use 
params.m.omega=.96;                                   %Source:
%% Tick Parameters
params.l.alpha=11.98;                                 %Source:
params.n.alpha=3.07;                                  %Source:
params.a.alpha=3.21;                                  %Source:
params.l.lambda=params.n.popsize*exp((params.l.alpha+3*params.n.alpha)/4); %Source:

%Beta's
% j=.0005;
% [c_l,c_n]=ContactRateCalc(params.m.popsize);
% % c_n=.002*params.m.popsize^.515*365;                             %Source:
% % c_l=.0013*params.m.popsize^.515*365;                               %Source:
% c_m=c_n*(params.n.popsize/params.m.popsize);                            %Source:
% params.n.beta=.967957*c_n*j;
% params.l.beta=.967957*c_l*j;
% params.m.beta=.5*c_m*j;
params.n.beta=2.0353;                                              %Derived from Beta estimations where infected tick proportion was ~.252 (Vaccine field data)
params.l.beta=1.3230;
params.m.beta=21.0270;

%% Time Params
params.tmax=30000; params.to=0;          %Total time the program runs
params.t=params.to:params.tmax;          %Set as 1/4 so that you can caluculate populations at each season for non-nymphal ticks
params.tol=.0000001;                     %Sets a tolerance that the model can run until instead of a set amount of time
                                         %if you wish tmax to be the only
                                         %qualifier set tol=0
                                         %In Solution function tmax is
                                         %hardcoded to be required so as to
                                         %prevent infinite loops for
                                         %unstable equilibrium 

end

