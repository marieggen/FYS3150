clear all
close all
clc

n = linspace(1,100,100);
%t = linspace(0.001,0.3,10);
t = 0.5;
x = linspace(0,1,100);
u = 0;

for j=1:length(t)
    for i=1:length(n)
        A = -(2.0/pi)*(1.0/n(i));
        term = A.*sin(n(i).*pi.*x).*exp(-((n(i)*pi)^2)*t(j));
        u = u + term;
    end
    plot(x,u+(1-x))
    hold('on')
    u = 0;
end



