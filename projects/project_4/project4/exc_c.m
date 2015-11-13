clear all
close all
clc

file_mcc1 = 'num_mcc_T1p0_Sallup.txt';%Temp = 1.0
file_mcc2 = 'num_mcc_T2p4_Srand.txt';%Temp = 2.4
num1 = load(file_mcc1(:));
num2 = load(file_mcc2(:));


numE1 = num1(:,1);
numCv1 = num1(:,2);
numAbsM1 = num1(:,5);
numAbsX1 = num1(:,6);
MCC1 = num1(:,8);
count1 = num1(:,9);

% figure(1)
% plot(MCC1,numE1,'k')
% title('Mean energy for T=1.0','fontsize',18)
% xlabel('MC cycles','fontsize',18)
% ylabel('<E/J>','fontsize',18)
% set(gca,'FontSize',15)
% 

% 
% figure(3)
% plot(MCC1,numAbsM1,'k')
% title('Mean absolute magnetization for T=1.0','fontsize',18)
% xlabel('MC cycles','fontsize',18)
% ylabel('<|M|>','fontsize',18)
% set(gca,'FontSize',15)


numE2 = num2(:,1);
numCv2 = num2(:,2);
numAbsM2 = num2(:,5);
numAbsX2 = num2(:,6);
MCC2 = num2(:,8);
count2 = num2(:,9)

% figure(4)
% plot(MCC2,numE2,'k')
% title('Mean energy for T=2.4','fontsize',18)
% xlabel('MC cycles','fontsize',18)
% ylabel('<E/J>','fontsize',18)
% set(gca,'FontSize',15)
% 

% 
% figure(6)
% plot(MCC2,numAbsM2,'k')
% title('Mean absolute magnetization for T=2.4','fontsize',18)
% xlabel('MC cycles','fontsize',18)
% ylabel('<|M|>','fontsize',18)
% set(gca,'FontSize',15)


len = length(count1);
M = 1.0;
for i=2:len
    count1(i) = count1(i)/M;
    count2(i) = count2(i)/M;
    M = M + 5.0;
end

mean(count1)
mean(count2)

figure(2)
plot(MCC1,count1,'k')
title('Accepted configurations for T=1.0','fontsize',18)
xlabel('MC cycles','fontsize',18)
ylabel('counts','fontsize',18)
set(gca,'FontSize',15)


figure(5)
plot(MCC2,count2,'k')
title('Accepted configurations for T=2.4','fontsize',18)
xlabel('MC cycles','fontsize',18)
ylabel('counts','fontsize',18)
set(gca,'FontSize',15)









