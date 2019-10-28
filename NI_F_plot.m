clear;
params=Lymeparams();
n_popsize=params.n.popsize;
m_popsize=params.m.popsize;
psi=0;
omega=params.m.omega;
mu=params.m.mu;
m_lambda=params.m.lambda;
BM=[11.85, 7.05, 1.00];
BL=[5.73, 3.41, .5];
ni=0:700;
F(4,1:701)=0;
figure
for k=1:701
    NI=ni(k);
F(1,k)=log(1-NI/n_popsize)+(BL(1)*m_popsize*exp(-psi*omega/4)*exp(-mu/4))/(4*exp(-mu/4)*m_popsize+m_lambda/4)*(1-exp(-BM(1)/2*NI/n_popsize))/(1-exp(-mu-psi*omega/4-BM(1)/2*NI/n_popsize));
F(2,k)=log(1-NI/n_popsize)+(BL(1)*m_popsize*exp(-psi*omega/4)*exp(-mu/4))/(4*exp(-mu/4)*m_popsize+m_lambda/4)*(1-exp(-BM(2)/2*NI/n_popsize))/(1-exp(-mu-psi*omega/4-BM(2)/2*NI/n_popsize));
F(3,k)=log(1-NI/n_popsize)+(BL(1)*m_popsize*exp(-psi*omega/4)*exp(-mu/4))/(4*exp(-mu/4)*m_popsize+m_lambda/4)*(1-exp(-BM(3)/2*NI/n_popsize))/(1-exp(-mu-psi*omega/4-BM(3)/2*NI/n_popsize));
end
plot(ni,F);
legend('\beta_M=11.85 \beta_L=5.73','\beta_M=7.05 \beta_L=3.41','\beta_M=1.00 \beta_L=.50')
ylabel('G(N_I)')
xlabel('N_I')