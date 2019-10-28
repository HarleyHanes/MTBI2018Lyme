function [sol,R] = LymeSolutions(params)
%LymeSolutions takes parameters and intial conditions from Lymeparams and
%outputs vectors of populations with respect to time
%   Detailed explanation goes here
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
to=params.to;
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
            elseif sum(((sol(6:10,end)-sol(6:10,end-1)))>tol)~=0
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
            kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
            a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
            b=exp(-params.m.mu);
            R=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
            j;
end

% Junk Code
% sol(6,:)=sol(1,:)/(sol(1,:)+sol(2,:)+sol(3,:));
% sol(7,:)=sol(2,:)/(sol(1,:)+sol(2,:)+sol(3,:));
% sol(8,:)=sol(3,:)/(sol(1,:)+sol(2,:)+sol(3,:));
% sol(9,:)=sol(4,:)/(sol(4,:)+sol(5,:));
% sol(10,:)=sol(5,:)/(sol(4,:)+sol(5,:));





%for k=0:20*100
%    psi=k/100*365;
%     for i=2:(tmax+1-to)
% %         if i==2
% %                 sol(1,i)=sol(1,i-1)*exp(-mu)*exp(-psi*omega/4)*exp((-beta_m/2*sol(5,i-1))/(sol(4,i-1)+sol(5,i-1)))+lambda_m/2*(exp(-3*mu/4)+3);
% %                 sol(2,i)=sol(2,i-1)*exp(-mu)+sol(1,i-1)*exp(-mu)*exp(-psi*omega/4)*(1-exp(-beta_m/2*sol(5,i-1)/(sol(4,i-1)+sol(5,i-1))));
% %                 sol(3,i)=sol(3,i-1)*exp(-mu)+sol(1,i-1)*exp(-mu)*(1-exp(-psi*omega/4));
% %                 sol(4,i)=lambda_l*exp(-(alpha_l+3*alpha_n)/4)*exp((-beta_l/4)*exp(-3*mu/4)*sol(2,i)/(sol(1,i)+sol(2,i)+sol(3,i)));
% %                 sol(5,i)=lambda_l*exp(-(alpha_l+3*alpha_n)/4)-sol(4,i);
% %                 sol(6,i)=sol(1,i)/(sol(1,i)+sol(2,i)+sol(3,i));             %6th row of sol plots % mice susceptible
% %                 sol(7,i)=sol(2,i)/(sol(1,i)+sol(2,i)+sol(3,i));             %7th row of sol plots % mice infected
% %                 sol(8,i)=sol(3,i)/(sol(1,i)+sol(2,i)+sol(3,i));             %8th row of sol plots % mice vaccinated
% %                 sol(9,i)=sol(4,i)/(sol(4,i)+sol(5,i));                      %8th row of sol plots % nymphs susceptible
% %                 sol(10,i)=sol(5,i)/(sol(4,i)+sol(5,i));                     %9th row of sol plots % nymphs infected
% %         else
% %                while sum(((sol(6:10,end)-sol(6:10,end-1))/sol(6:10,end-1))>tol)~=0
%                     sol(1,i)=sol(1,i-1)*exp(-mu)*exp(-psi*omega/4)*exp((-beta_m/2*sol(5,i-1))...
%                         /(sol(4,i-1)+sol(5,i-1)))+lambda_m/4*(exp(-3*mu/4)+3);
%                     sol(2,i)=sol(2,i-1)*exp(-mu)+sol(1,i-1)*exp(-mu)*exp(-psi*omega/4)*...
%                         (1-exp(-beta_m/2*sol(5,i-1)/(sol(4,i-1)+sol(5,i-1))));
%                     sol(3,i)=sol(3,i-1)*exp(-mu)+sol(1,i-1)*exp(-mu)*(1-exp(-psi*omega/4));
%                     sol(4,i)=lambda_l*exp(-(alpha_l+3*alpha_n)/4)*exp(((-beta_l/4)*exp(-mu/4)*...
%                         sol(2,i-1)+sol(1,i-1)*exp(-mu/4)*exp(-psi*omega/4)*...
%                         (1-exp(-(beta_m/2)*sol(5,i-1)/(sol(5,i-1)+sol(4,i-1)))))/...
%                         (exp(-mu/4)*(sol(1,i-1)+sol(2,i-1)+sol(3,i-1))+lambda_m/4));
%                     sol(5,i)=lambda_l*exp(-(alpha_l+3*alpha_n)/4)-sol(4,i);
%                     sol(6,i)=sol(1,i)/(sol(1,i)+sol(2,i)+sol(3,i));             %6th row of sol plots % mice susceptible
%                     sol(7,i)=sol(2,i)/(sol(1,i)+sol(2,i)+sol(3,i));             %7th row of sol plots % mice infected
%                     sol(8,i)=sol(3,i)/(sol(1,i)+sol(2,i)+sol(3,i));             %8th row of sol plots % mice vaccinated
%                     sol(9,i)=sol(4,i)/(sol(4,i)+sol(5,i));                      %8th row of sol plots % nymphs susceptible
%                     sol(10,i)=sol(5,i)/(sol(4,i)+sol(5,i));                     %9th row of sol plots % nymphs infected
% %                end
%         end
%          eq(1,k)=sol(1,end);
%          eq(2,k)=sol(2,end);
%          eq(3,k)=sol(3,end);
%          eq(4,k)=sol(4,end);
%          eq(5,k)=sol(5,end);
%          eq(6,k)=sol(6,end);
%          eq(7,k)=sol(7,end);
%          eq(8,k)=sol(8,end);
%          eq(9,k)=sol(9,end);
%          eq(10,k)=sol(10,end);
%     end


