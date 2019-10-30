function [sol,R] = LymeSolutions(varargin)
%LymeSolutions takes parameters and intial conditions from Lymeparams and
%outputs vectors of populations with respect to time
%   Detailed explanation goes here
%Load in inputs
if nargin==1
    params=varargin{1};
    toltest=1;
elseif nargin==2
    params=varargin{1};
    toltest=varagin{2};          %if ==0 then won't end solver at tol.
else
    warning("LymeSolutions only takes 1 or 2 inputs, %i inputs given",nargin)
    keyboard
end
%% Redefining params for concise coding, see Lymeparams for definitions
%Mouse
beta_m=params.m.beta;
mu=params.m.mu;
lambda_m=params.m.lambda;
psi=params.m.psi;
omega=params.m.omega;
%Tick
alpha_l=params.l.alpha;
lambda_l=params.l.lambda;
beta_l=params.l.beta;
alpha_n=params.n.alpha;
beta_n=params.n.beta;
%a.alpha=params.a.alpha;
%% Intial Conditions
% sol(6,1)=sol(1,1)/(sol(1,1)+sol(2,1)+sol(3,1));             %6th row of sol plots % mice susceptible
% sol(7,1)=sol(2,1)/(sol(1,1)+sol(2,1)+sol(3,1));             %7th row of sol plots % mice infected
% sol(8,1)=sol(3,1)/(sol(1,1)+sol(2,1)+sol(3,1));             %8th row of sol plots % mice vaccinated
% sol(9,1)=sol(4,1)/(sol(4,1)+sol(5,1));                      %8th row of sol plots % nymphs susceptible
% sol(10,1)=sol(5,1)/(sol(4,1)+sol(5,1));                     %9th row of sol plots % nymphs infected

%Time and tolerance
tmax=params.tmax;
%to=params.to;
tol=params.tol;
sol(1,1)=params.m.So;           %First row of sol plots susceptible mice
sol(2,1)=params.m.Io;           %Second row of sol plots infected mice
sol(3,1)=params.m.Vo;           %Third row of sol plots vaccinated mice
sol(4,1)=params.n.So;           %Fourth row of sol plots susceptible nymphs
sol(5,1)=params.n.Io;           %Fifth row of sol plots infected nymphs
    for i=1:tmax
            if i==1
                %Defines Modular Probabilities
                    p_m_noinfect=exp(-beta_m/2*sol(5,i)/(sol(4,i)+sol(5,i)));
                    p_m_novacc=exp(-psi*omega/4);
                    p_m_nodeath=exp(-mu);
                    m_recruit=lambda_m/4*(exp(-3*mu/4)+3);
                    n_recruit=lambda_l*exp(-(alpha_l+3*alpha_n)/4);
                        n_prop_miceinfect=(sol(2,i)+sol(1,i)*p_m_novacc*(1-p_m_noinfect))*exp(-mu/4);
                        n_prop_micetotal=exp(-mu/4)*(sol(1,i)+sol(2,i)+sol(3,i))+lambda_m/4;
                    p_n_noinfect=exp((-beta_l/4)*n_prop_miceinfect/n_prop_micetotal);
            %Population Totals
            sol(1,i+1)=sol(1,i)*p_m_nodeath*p_m_novacc*p_m_noinfect+m_recruit;
            sol(2,i+1)=sol(2,i)*p_m_nodeath+sol(1,i)*p_m_nodeath*p_m_novacc*(1-p_m_noinfect);
            sol(3,i+1)=sol(3,i)*p_m_nodeath+sol(1,i)*p_m_nodeath*(1-p_m_novacc);
            sol(4,i+1)=n_recruit*p_n_noinfect;
            sol(5,i+1)=n_recruit*(1-p_n_noinfect); 
            sol(6,i+1)=sol(1,i+1)/(sol(1,i+1)+sol(2,i+1)+sol(3,i+1));             %6th row of sol plots % mice susceptible
            sol(7,i+1)=sol(2,i+1)/(sol(1,i+1)+sol(2,i+1)+sol(3,i+1));              %7th row of sol plots % mice infected
            sol(8,i+1)=sol(3,i+1)/(sol(1,i+1)+sol(2,i+1)+sol(3,i+1));             %8th row of sol plots % mice vaccinated
            sol(9,i+1)=sol(4,i+1)/(sol(4,i+1)+sol(5,i+1));                        %8th row of sol plots % nymphs susceptible
            sol(10,i+1)=sol(5,i+1)/(sol(4,i+1)+sol(5,i+1));    %9th row of sol plots % nymphs infected  
            j=1;
            elseif norm(sol(6:10,end)-sol(6:10,end-1))>tol || toltest==0
            %Defines Modular Probabilities
                    p_m_noinfect=exp(-beta_m/2*sol(5,i)/(sol(4,i)+sol(5,i)));
                    p_m_novacc=exp(-psi*omega/4);
                    p_m_nodeath=exp(-mu);
                    m_recruit=lambda_m/4*(exp(-3*mu/4)+3);
                    n_recruit=lambda_l*exp(-(alpha_l+3*alpha_n)/4);
                        n_prop_miceinfect=(sol(2,i)+sol(1,i)*p_m_novacc*(1-p_m_noinfect))*exp(-mu/4);
                        n_prop_micetotal=exp(-mu/4)*(sol(1,i)+sol(2,i)+sol(3,i))+lambda_m/4;
                    p_n_noinfect=exp((-beta_l/4)*n_prop_miceinfect/n_prop_micetotal);
            %Population Totals
            sol(1,i+1)=sol(1,i)*p_m_nodeath*p_m_novacc*p_m_noinfect+m_recruit;
            sol(2,i+1)=sol(2,i)*p_m_nodeath+sol(1,i)*p_m_nodeath*p_m_novacc*(1-p_m_noinfect);
            sol(3,i+1)=sol(3,i)*p_m_nodeath+sol(1,i)*p_m_nodeath*(1-p_m_novacc);
            sol(4,i+1)=n_recruit*p_n_noinfect;
            sol(5,i+1)=n_recruit*(1-p_n_noinfect); 
            sol(6,i+1)=sol(1,i+1)/(sol(1,i+1)+sol(2,i+1)+sol(3,i+1));             %6th row of sol plots % mice susceptible
            sol(7,i+1)=sol(2,i+1)/(sol(1,i+1)+sol(2,i+1)+sol(3,i+1));              %7th row of sol plots % mice infected
            sol(8,i+1)=sol(3,i+1)/(sol(1,i+1)+sol(2,i+1)+sol(3,i+1));             %8th row of sol plots % mice vaccinated
            sol(9,i+1)=sol(4,i+1)/(sol(4,i+1)+sol(5,i+1));                        %8th row of sol plots % nymphs susceptible
            sol(10,i+1)=sol(5,i+1)/(sol(4,i+1)+sol(5,i+1));  %9th row of sol plots % nymphs infected
            j=j+1;
            end
    end
    R=CalculateR(params);
end
