function [c_larva,c_nymph] = ContactRateCalc(m_popsize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

TempData=xlsread('temp_dutchess_2008_2018.xls','E2208:F3302');
%TempData=tempdutchess3;                 %Import Rows 2208-3302, Temp Min and Max Columns
TempA=(TempData(:,1)+TempData(:,2))/2;
t=1:length(TempData(:,1));
%plot(t,TempData(:,1))
%plot(t,TempData(:,2))
a_1=.1057; b_1 = 0.0374; a_2 = 0.0481; b_2 = 0.04018; a_3 = 0.009363; b_3 = 0.01734;
a_4 = 0.001915; b_4 = 0.001805; a_5 = 0.004033; b_5 = 0.003589; a_6 = 0.001028;
b_6 = 0.006714; a_7 = 0.001692; b_7 = .002853; c = 0.06712;

a=[a_1 a_2 a_3 a_4 a_5 a_6 a_7];
b=[b_1 b_2 b_3 b_4 b_5 b_6 b_7];
n=1:7;
for i=1:length(TempA)
theta_t(i)=c+sum(a.*cos(2*n.*pi.*TempA(i)/365)+b.*sin(2.*n.*pi.*TempA(i)/365));
end

for i=1:3
    for k=120:273
        thetaA(i)=mean(theta_t((i-1)*365+k)*360);
    end
end
theta_mean=mean(thetaA);
Mstar=m_popsize*log(2/.0012);
c_larva=.0013*Mstar^.485*theta_mean*365;
c_nymph=.002*Mstar^.485*theta_mean*365;
end

