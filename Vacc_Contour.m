function [Beta_c,Psi] = Vacc_Contour(params, contour)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
c=40;                                         %Number of iterations on Psi
Psi=linspace(0,.06,c);
Beta_c=zeros(length(contour),c);
for i=1:c
    params.m.psi=Psi(i)*365;
    [eq,Beta]=BetaEstimation(params);
    B=Beta(1,:);
    for j=1:length(contour)
        percent=contour(j);
        if sum(eq(5,:)>=percent)~=0
            Beta_c(j,i)=min(B(eq(5,:)>=percent));
        end
    end
end

        

end

