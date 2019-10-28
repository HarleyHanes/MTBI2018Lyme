clear; close all
%% Lyme Disease Main script file
% Run all main programs through this script file to produce required plots
%Sections
%   Section 1: Population vs. Time
%   Section 2: Beta Estimation
%Required Files
%   Lymeparams.m            
%   VaccineEffect.m
%   VaccineEffect_Single.m
%   Vacc_Contour.m
%   PsiEstimation.m
%   LymeSolutions.m
%   Cost_Analysise
%   Cost_Analysis_Year.m
%   Cost_Analysis_Psi.m
%   ContactRateCalc.m
%   BetaEstimation.m
%   BetaCalc.m


%% Section 1: Population vs. Time
% This section of code creates a plot of the mice and tick populations with
% respected to time at parameter values selected from Lymeparams.
% Subfiles: Lymeparams.m, Lyme Solutions.m

%Parameter Loading
%Use Lymeparams to load correct parameter values
params=Lymeparams();
%Modify Infection rates
params.n.beta=1.47;
params.l.beta=params.n.beta/.2;
params.m.beta=params.n.beta/.0968;
%Set vaccination or not
params.m.psi=10;  %With vaccination
%params.m.psi=0;   %Without vaccination

%Use LymeSolutions to create matrix of compartment sizes vs. time and
%create time vector
[sol,R]=LymeSolutions(params);
t=0:(length(sol(1,:))-1);

%Create Plots
figure
%Mouse Population Plots
subplot(2,1,1)
plot(t,sol(1,:),'^','MarkerFaceColor','b')
hold on
plot(t,sol(2,:),'s','MarkerFaceColor','r')
plot(t,sol(3,:),'*','MarkerFaceColor','y')
%plot(t,sol(3,:),'s')
hold off
legend({'Susceptible', 'Infected','Vaccinated'},'Fontsize',12)
%ylabel({'Number','of Mice'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 12)
%title('Mouse Population','FontSize')
axis([0 6 0 50])

%Tick Population Plots
subplot(2,1,2)
plot(t,sol(4,:),'^','MarkerFaceColor','b')
hold on
plot(t,sol(5,:),'s','MarkerFaceColor','r')
hold off
%legend('Susceptible Ticks','Infected Ticks')
%ylabel({'Number','of Ticks'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
xlabel('Time (years)','FontSize',12)

%% WHAT IS THIS???
kappa=exp((-params.l.alpha-3*params.n.alpha)/4)*(1-exp(-params.m.mu))*params.l.beta*params.l.lambda/(3*exp(-params.m.mu)+exp(-3*params.m.mu/4)*params.m.lambda);
a=params.m.beta*params.m.lambda/(8*params.l.lambda)*((3*exp(-params.m.mu/4)+exp(-params.m.mu))*exp(-3*params.m.mu/4-params.m.psi*params.m.omega/4))/((1-exp(-params.m.mu-params.m.psi*params.m.omega/4))*exp((-params.l.alpha-3*params.n.alpha)/4));
b=exp(-params.m.mu);
R=.5*(kappa*a+sqrt((kappa*a)^2+4*kappa*a*b/(1-b)));
xt = get(gca, 'XTick');
set(gca, 'FontSize', 12);
axis([0 6 0 1000])
% figure
% 
% plot(t,sol(6,:))
% hold on
% plot(t,sol(7,:))
% plot(t,sol(8,:))
% plot(t,sol(9,:))
% plot(t,sol(10,:))
% hold off
% legend('% mice susceptible','% mice infectious','% mice vaccinated','% tick susceptible','% tick infected')
% 
% r=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu))

%% Section 2: Beta and R Estimation
% This section of code creates a plot of the equilibrium infection rates in 
% mice and ticks with changing Betas (x-axis is Beta_L see BetaEstimation.m 
% for explanation)
% Subfiles: Lymeparams.m, BetaEstimation.m

%Step 1: Reload paramters
params=Lymeparams();

%Step 2: Run BetaEstimation To get vectors of Beta, Equilibrium Values, and R
[eq,R,Beta]=BetaEstimation(params,1.4,200);

% Step 3: Form Plots
figure
plot(Beta(2,:),eq(7,:),'LineWidth',2)
grid on
hold on 
plot(Beta(2,:),eq(10,:),':','LineWidth',2)
%Mark 3 Beta values showing R0 at different equilibrium points
TestBetas=[Beta(2,98),Beta(2,123),Beta(2,164)];
TestEq=[eq(10,98),eq(10,123),eq(10,164)];
plot(TestBetas(1),TestEq(1),'d','MarkerFaceColor','k');
plot(TestBetas(2),TestEq(2),'^','MarkerFaceColor','k');
plot(TestBetas(3),TestEq(3),'s','MarkerFaceColor','k');
hold off
%Graph formating
legend({'Infectious Mice','Infectious Nymphs','R_0=3.02','R_0=4.77','R_0=8.51'},'Location','northwest','FontSize',12)
xlabel('\beta_L (Years^{-1})','FontSize',12)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 12);
ylabel({'Proportion of','Population'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
str=sprintf('Asymptotic Equilibria of Variable \\beta''s at, M=%.0f, N=%.0f, \\psi=%.0f',...
    params.m.popsize, params.n.popsize, params.m.psi);
title({str});
%figure
%plot(Beta(2,:),R)

%axis([0 .02 0 1])

%Plot equilibrium values vs. R
% plot(R,eq(7,:))
% hold on
% plot(R,eq(10,:))
% hold off
% legend('Proportion of Mice Infectious','Proportion Nymphs Infectious')
% xlabel('r Stability Condition')
% ylabel('% of Total Population')
% %axis([1.3 1.4 0 .005])



%r=(params.m.beta*params.l.beta)/(8*(1-exp(-params.m.mu-params.m.psi*params.m.omega/4)))*(.25*exp(-params.m.mu)+.75*exp(-params.m.mu/4))/(.25*exp(-3*params.m.mu/4)+.75*exp(-params.m.mu))

%% Psi Estimation
clear;
%Step 1: Load Parameters
params=Lymeparams();

%Step 2: Run PsiEstimation
[eq,R,Psi]=PsiEstimation(params);

%Step 3: Plot Results
figure
i=1;
R1=R(1,:);
R2=R(2,:);
R3=R(3,:);
P1=min(Psi(R(1,:)<1));
P2=min(Psi(R(2,:)<1));
P3=min(Psi(R(3,:)<1));
semilogy(Psi,R(1,:),'LineWidth',2)
grid on
hold on
semilogy(Psi,R(2,:),':','LineWidth',2)
semilogy(Psi,R(3,:),'--','LineWidth',2)
semilogy(P1,1,'d','MarkerFaceColor','k')
semilogy(P2,1,'^','MarkerFaceColor','k')
semilogy(P3,1,'s','MarkerFaceColor','k')
hold off
xlabel('\psi (Years^{-1})','FontSize',12)
ylabel('R_c','Rotation',0,'FontSize',12,'HorizontalAlignment','right')
%title('Disease Free Equilibrium Stability as a Function of Psi')
legend({'\beta_N=1.47, \beta_L=5.73, \beta_M=11.85','\beta_N=.86, \beta_L=4.29, \beta_M=8.87','\beta_N=.68, \beta_L=3.41,\beta_M=7.05','\psi=4.58','\psi=6.55','\psi=11.04'},'FontSize',12)
axis([0 12 -.00001 50])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12);

for i=1
figure
semilogy(Psi,eq.S(7,:),'LineWidth',2)
grid on
hold on
semilogy(Psi,eq.S(8,:),'--','LineWidth',2)
semilogy(Psi,eq.S(10,:),':','LineWidth',2)
hold off
%title({'Asymptotic Equilibria at tol=.0000001, M=50, N=1000','\beta_N=3.618, \beta_L=2.35,\beta_M=37.38'})
legend({'Infectious Mice','Vaccinated Mice','Infectious Nymphs'},'FontSize',12)
xlabel('\psi (Years^{-1})','FontSize',12)
ylabel({'Proportion of', 'Population'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
axis([0 10 .000001 1])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12);

figure
semilogy(Psi,eq.M(7,:),'LineWidth',2)
grid on
hold on
semilogy(Psi,eq.M(8,:),'--','LineWidth',2)
semilogy(Psi,eq.M(10,:),':','LineWidth',2)
hold off
%title({'Asymptotic Equilibria at tol=.0000001, M=50, N=1000','\beta_N=17.19, \beta_L=11.17,\beta_M=117.54'})
legend({'Infectious Mice','Vaccinated Mice','Infectious Nymphs'},'FontSize',12)
xlabel('\psi (Years^{-1})','FontSize',12)
ylabel({'Proportion of','Population'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
axis([0 10 .000001 1])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12);

figure
semilogy(Psi,eq.B(7,:),'LineWidth',2)
grid on
hold on
semilogy(Psi,eq.B(8,:),'--','LineWidth',2)
semilogy(Psi,eq.B(10,:),':','LineWidth',2)
hold off
%title({'Asymptotic Equilibria at tol=.0000001, M=50, N=1000','\beta_N=36.63, \beta_L=23.81,\beta_M=378.44'})
legend({'Infectious Mice','Vaccinated Mice','Infectious Nymphs'},'FontSize',12)
xlabel('\psi (Years^{-1})','FontSize',12)
ylabel({'Proportion of','Population'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
axis([0 10 .000001 1])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12);
 end

%% Vaccine Effectiveness
clear
year=[2 5 10];
params=Lymeparams();

[eq,R,Psi]=VaccineEffect(params,year);
%% Plots
for i=1
%New Plots
figure
subplot(1,3,1)
plot(Psi,eq.S.nymph(1,:),'LineWidth',2)
hold on
plot(Psi,eq.S.nymph(2,:),'--','LineWidth',2)
plot(Psi,eq.S.nymph(3,:),':','LineWidth',2)
plot(Psi,eq.S.nymph(4,:),'-.','LineWidth',2)
axis([0 12 0 .6])
title('\beta_N=.68, \beta_L=3.41,\beta_M=7.05','Fontsize',12)
ylabel({'Proportion','of Nymph','Population'},'Rotation',0,'FontSize',12,'HorizontalAlignment','right')
legend({'No Vaccination','2 years','5 years', '10 years'},'FontSize',12)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12);
    
subplot(1,3,2)
plot(Psi,eq.M.nymph(1,:),'LineWidth',2)
hold on
plot(Psi,eq.M.nymph(2,:),'--','LineWidth',2)
plot(Psi,eq.M.nymph(3,:),':','LineWidth',2)
plot(Psi,eq.M.nymph(4,:),'-.','LineWidth',2)
axis([0 12 0 .6])
title('\beta_N=.86, \beta_L=4.29,\beta_M=8.87','Fontsize',12)
xlabel('\psi (Years^{-1})','FontSize',12)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 12);
subplot(1,3,3)
plot(Psi,eq.B.nymph(1,:),'LineWidth',2)
hold on
plot(Psi,eq.B.nymph(2,:),'--','LineWidth',2)
plot(Psi,eq.B.nymph(3,:),':','LineWidth',2)
plot(Psi,eq.B.nymph(4,:),'-.','LineWidth',2)
axis([0 12 0 .6])
title('\beta_N=1.47, \beta_L=5.73,\beta_M=11.85','Fontsize',12)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 12);
% legend('No Vaccination','2 years','5 years', '10 years')
end

% Old Plots
for i=1
% figure
% subplot(1,3,1)
% plot(Psi,eq.S.nymph(1,:))
% hold on
% plot(Psi,eq.S.nymph(2,:))
% plot(Psi,eq.S.nymph(3,:))
% plot(Psi,eq.S.nymph(4,:))
% axis([0 .06 0 1])
% title('% Nymphs Infected')
% ylabel('% of Total Population')
% 
% subplot(1,3,2)
% plot(Psi,eq.S.mouse(1,:))
% hold on
% plot(Psi,eq.S.mouse(2,:))
% plot(Psi,eq.S.mouse(3,:))
% plot(Psi,eq.S.mouse(4,:))
% axis([0 .06 0 1])
% xlabel('\psi (days^{-1})')
% title({'Vaccine Effectiveness Compared to Years of Use, M=50, N=1000 \beta_N=.0028, \beta_L=.0018,\beta_M=.0286','% Mice Infected'})
% 
% subplot(1,3,3)
% plot(Psi,eq.S.vacc(1,:))
% hold on
% plot(Psi,eq.S.vacc(2,:))
% plot(Psi,eq.S.vacc(3,:))
% plot(Psi,eq.S.vacc(4,:))
% hold off
% legend('No Vaccination','2 years','5 years', '10 years')
% title('% Mice Vaccinated')
% axis([0 .06 0 1])
% hold off
% 
% figure
% subplot(1,3,1)
% plot(Psi,eq.M.nymph(1,:))
% hold on
% plot(Psi,eq.M.nymph(2,:))
% plot(Psi,eq.M.nymph(3,:))
% plot(Psi,eq.M.nymph(4,:))
% axis([0 .06 0 1])
% title('% Nymphs Infected');ylabel('% of Total Population')
% 
% subplot(1,3,2)
% plot(Psi,eq.M.mouse(1,:))
% hold on
% plot(Psi,eq.M.mouse(2,:))
% plot(Psi,eq.M.mouse(3,:))
% plot(Psi,eq.M.mouse(4,:))
% axis([0 .06 0 1])
% xlabel('\psi (days^{-1})')
% title({'Vaccine Effectiveness Compared to Years of Use, M=50, N=1000 \beta_N=.0417, \beta_L=.0271,\beta_M=.4304','% Mice Infected'})
% 
% subplot(1,3,3)
% plot(Psi,eq.M.vacc(1,:))
% hold on
% plot(Psi,eq.M.vacc(2,:))
% plot(Psi,eq.M.vacc(3,:))
% plot(Psi,eq.M.vacc(4,:))
% hold off
% legend('No Vaccination','2 years','5 years', '10 years')
% title('% Mice Vaccinated')
% axis([0 .06 0 1])
% hold off
% 
% figure
% subplot(1,3,1)
% plot(Psi,eq.B.nymph(1,:))
% hold on
% plot(Psi,eq.B.nymph(2,:))
% plot(Psi,eq.B.nymph(3,:))
% plot(Psi,eq.B.nymph(4,:))
% axis([0 .06 0 1])
% title('% Nymphs Infected')
% ylabel('% of Total Population')
% 
% subplot(1,3,2)
% plot(Psi,eq.B.mouse(1,:))
% hold on
% plot(Psi,eq.B.mouse(2,:))
% plot(Psi,eq.B.mouse(3,:))
% plot(Psi,eq.B.mouse(4,:))
% axis([0 .06 0 1])
% xlabel('\psi (days^{-1})')
% title({'Vaccine Effectiveness Compared to Years of Use, M=50, N=1000 \beta_N=.0833, \beta_L=.0542,\beta_M=.8609','% Mice Infected'})
% 
% subplot(1,3,3)
% plot(Psi,eq.B.vacc(1,:))
% hold on
% plot(Psi,eq.B.vacc(2,:))
% plot(Psi,eq.B.vacc(3,:))
% plot(Psi,eq.B.vacc(4,:))
% hold off
% legend('No Vaccination','2 years','5 years', '10 years')
% title('% Mice Vaccinated')
% axis([0 .06 0 1])
% hold off
end

%% Contour plots
contour=[100, 500, 700];
params=Lymeparams();
[Beta_c, Psi]=Vacc_Contour(params,contour);
Beta_c=Beta_c/365;
%%
plot(Psi,Beta_c)
legend
%axis([0 .06 0 500])

%% Cost Analysis-Psi
clear;
H_pop=[80 160 750];
year=10;
for kk=1:length(H_pop)
    figure
    grid on
    gamma=[.1; .25; .5; 1];
    H_popsize=H_pop(kk);
    for i=1:length(gamma)  
        params=Lymeparams();
        [C_v,C_c,DALY,Psi,H]=Cost_Analysis_Psi(params,H_popsize,gamma(i),year);
        C_vacc(i,:)=C_v;
        C_case(i,:)=C_c;
        DALY_time(i,:)=DALY;
        Humans(i,:)=H;
    end
    C_total=C_vacc+C_case;
    C_change=C_total(:,1)-C_total;
    DALY_saved=DALY_time(:,1)-DALY_time;
    for i=1:length(gamma)
        if sum(C_change(i,:)<0)>0
            C_zeros=C_change(C_change(i,:)<0);
            d=length(C_change(i,:))-length(C_zeros);
            Equals(i)=Psi(d);
        end
    end                     % Finding points where you lose money
    for i=1:length(gamma)
        delta=C_change(i,2:end)-C_change(i,1:end-1);
        d=length(delta)-length(delta(delta<0))+1;
        Inflect(1,i)=Psi(d);
        Inflect(2,i)=C_change(i,d);
    end
    hold on
    %legend('\gamma=.183','\gamma=.456','\gamma=.913','\gamma=1.825')
    if length(Inflect(Inflect(1,:)>0))==1
        plot(Inflect(1,4),Inflect(2,4),'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,4)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,4)))])
    elseif length(Inflect(Inflect(1,:)>0))==2
        plot(Inflect(1,4),Inflect(2,4),'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,4)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,4)))])
        plot(Inflect(1,3),Inflect(2,3),'^','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,3)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,3)))])
    elseif length(Inflect(Inflect(1,:)>0))==3
        plot(Inflect(1,4),Inflect(2,4),'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,4)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,4)))])
        plot(Inflect(1,3),Inflect(2,3),'^','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,3)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,3)))])
        plot(Inflect(1,2),Inflect(2,2),'d','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,2)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,2)))])
    %    legend(['\psi= ', num2str(Equals(1))],['\psi= ', num2str(Equals(2))], ['\psi= ', num2str(Equals(3))])
    elseif length(Inflect(Inflect(1,:)>0))==4
        plot(Inflect(1,4),Inflect(2,4),'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,4)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,4)))])
        plot(Inflect(1,3),Inflect(2,3),'^','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,3)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,3)))])
        plot(Inflect(1,2),Inflect(2,2),'d','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,2)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,2)))])
        plot(Inflect(1,1),Inflect(2,1),'o','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(round(Inflect(1,1)*100)/100)])%,', Savings=$', num2str(round(Inflect(2,1)))])
    end
%    if kk==1
        plot(Psi,C_change(1,:),'k','DisplayName',['\gamma= ', num2str(gamma(1))],'LineWidth',2)
        plot(Psi,C_change(2,:),':b','DisplayName',['\gamma= ', num2str(gamma(2))],'LineWidth',2)
        plot(Psi,C_change(3,:),'--g','DisplayName',['\gamma= ', num2str(gamma(3))],'LineWidth',2)
        plot(Psi,C_change(4,:),'-.r','DisplayName',['\gamma= ', num2str(gamma(4))],'LineWidth',2)
%     else
%         plot(Psi,C_change(1,:),'k','LineWidth',2,'HandleVisibility','off')
%         plot(Psi,C_change(2,:),':b','LineWidth',2,'HandleVisibility','off')
%         plot(Psi,C_change(3,:),'--g','LineWidth',2,'HandleVisibility','off')
%         plot(Psi,C_change(4,:),'-.r','LineWidth',2,'HandleVisibility','off')
%     end
    hold off
    if kk==2
        title('High Traffic Trail','FontSize',15)
    elseif kk==1
        title('Low Traffic Trail','FontSize',15)
     elseif kk==3
        title('Suburban Park','FontSize',15)
    end
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 14);
    xlabel('\psi (Year^{-1})','FontSize',14)
    ylabel({'Dollars',' Saved','(Year^{-1})'},'FontSize',14,'Rotation',0,'HorizontalAlignment','right')
    legend({},'FontSize',15)
end
        % if length(Equals)==1
    %     plot(Equals(1),0,'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(1))])
    % elseif length(Equals)==2
    %     plot(Equals(1),0,'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(1))])
    %     plot(Equals(2),1,'^','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(2))])
    % elseif length(Equals)==3
    %     plot(Equals(1),0,'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(1))])
    %     plot(Equals(2),1,'^','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(2))])
    %     plot(Equals(3),1,'d','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(3))])
    % %    legend(['\psi= ', num2str(Equals(1))],['\psi= ', num2str(Equals(2))], ['\psi= ', num2str(Equals(3))])
    % elseif length(Equals)==4
    %     plot(Equals(1),0,'s','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(1))])
    %     plot(Equals(2),1,'^','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(2))])
    %     plot(Equals(3),1,'d','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(3))])
    %     plot(Equals(4),1,'o','MarkerFaceColor','k','DisplayName',['\psi= ', num2str(Equals(4))])
    % %    legend(['\psi= ', num2str(Equals(1))],['\psi= ', num2str(Equals(2))],['\psi= ', num2str(Equals(3))],['\psi= ', num2str(Equals(4))])
    % end
%yyaxis right
%plot(Psi,DALY_saved)
%legend('\gamma=.183','\gamma=.456','\gamma=.913','\gamma=1.825')
%title('M=50, N=1000, \beta_N=.0028, \beta_L=.0018,\beta_M=.0286'), axis([0 .006 0 7000])
%title('10 years, H=500, M=50, N=1000, \beta_N=.0417, \beta_L=.0271,\beta_M=.4304')
%title('10 years, H_s=200, M=50, N=1000, \beta_N=.0833, \beta_L=.0542,\beta_M=.8609')
%axis([0 .13*365/4 -20000 200000])
%%
%Effect=Humans;
for i=1:4
    per_person(i,:)=DALY_saved(i,:)/H_popsize(i)*year
end
plot(C_vacc(1,:),per_person(1,:),'LineWidth',2)
xlabel('Dollars Spent on Vaccination')
ylabel('DALYs Saved')
legend('2','5','10','15','20')
%axis([0 30000 0 60])
%axis([0 100000 0 60])

%% Cost Analysis-Years


