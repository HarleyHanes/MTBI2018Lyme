% t=0:.001:24*pi;
close all;
clear;
% x=sin(.25*t);
% y=-sin(.25*t);
% plot(t,x);
% hold on
% plot(t,y);
% hold off
% axis([0 24*pi 0 1.1])
t=0:.0001:60;
t2=-22.5:.0001:0;
x=normpdf(t2,-17.5,1);
for k=1:22.5/.0001+1
    if x(k)<.001
        x(k)=-2;
    end
end
plot(t2,x,'b--')
hold on
x1=normpdf(t,5,1);
for k=1:60/.0001+1
    if x1(k)<.001
        x1(k)=-2;
    end
end
plot(t,x1,'b--')

x2=normpdf(t,12.5,1);
for k=1:60/.0001+1
    if x2(k)<.001
        x2(k)=-2;
    end
end
plot(t,x2,'b')


x3=normpdf(t, 20,1);
for k=1:60/.0001+1
    if x3(k)<.001
        x3(k)=-2;
    end
end
plot(t,x3,'b--')

% x4=normpdf(t,27.5,1);
% plot(t,x4,'b')
x5=normpdf(t,35,1);
for k=1:60/.0001+1
    if x5(k)<.001
        x5(k)=-2;
    end
end
plot(t,x5,'b')

x6=normpdf(t,42.5,1);
for k=1:60/.0001+1
    if x6(k)<.001
        x6(k)=-2;
    end
end
plot(t,x6,'b--')
x7=normpdf(t,50,1);
for k=1:60/.0001+1
    if x7(k)<.001
        x7(k)=-2;
    end
end
plot(t,x7,'b')
hold off
% hold on
% 
% for i=1:8
%     x=normpdf(t,7.5*i,1);
%     plot(t,x)
% end
% hold off
% % plot(t,x1)
% % hold on
% % plot(t,x2)
% % plot(t,x3)
% % hold off
 axis([-22.5 60 0 .5])
 set(gca,'YTickLabel',[])
 set(gca,'XTickLabel',[])
