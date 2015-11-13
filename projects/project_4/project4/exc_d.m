clear all
close all 
clc

file = 'P_E_T2p4.txt';
P = load(file(:));
E = linspace(-800,800,1600);
plot(E,P,'r')
title('Probability of total energy','fontsize',18)
xlabel('E/J','fontsize',18)
ylabel('P(E/J)','fontsize',18)
set(gca,'FontSize',15)
hold('on')
x = linspace(-800,800,5000);
T = 2.4;
beta = -(1.0/T);
y = 4*normpdf(x,-490,54);
plot(x,y,'k','linewidth',3)
title('Probability of total energy','fontsize',18)
xlabel('E/J','fontsize',18)
ylabel('P(E/J)','fontsize',18)
set(gca,'FontSize',15)
legend('P(E)','Normal distribution','fontsize',18)


