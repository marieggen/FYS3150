clear all
close all
clc

% N=20, #MCC = 1e5
file_0 = 'num_20_0.txt';
file_1 = 'num_20_1.txt';
file_2 = 'num_20_2.txt';
file_3 = 'num_20_3.txt';

data = [load(file_0(:)) ; load(file_1(:)) ; load(file_2(:)) ; load(file_3(:))];

E = data(:,1);
Cv = data(:,2);
absM = data(:,5);
absX = data(:,6);
T = data(:,7);

% N=20, #MCC = 1e6
filemill_0_20 = 'millnum_20_0.txt';
filemill_1_20 = 'millnum_20_1.txt';
filemill_2_20 = 'millnum_20_2.txt';
filemill_3_20 = 'millnum_20_3.txt';

datamill1 = [load(filemill_0_20(:)) ; load(filemill_1_20(:)) ; load(filemill_2_20(:)) ; load(filemill_3_20(:))];

Emill20 = datamill1(:,1);
Cvmill20 = datamill1(:,2);
absMmill20 = datamill1(:,5);
absXmill20 = datamill1(:,6);

% figure(1)  
% plot(T,Cv,'ko-')
% hold('on')
% plot(T,Cvmill20,'ro-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<Cv/k>','fontsize', 18)
% title('<Cv/k> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% legend('#MCC = 1e5', '#MCC = 1e6','location','northwest')

% N=40, #MCC = 1e6
filemill_0_40 = 'millnum_40_0.txt';
filemill_1_40 = 'millnum_40_1.txt';
filemill_2_40 = 'millnum_40_2.txt';
filemill_3_40 = 'millnum_40_3.txt';

datamill2 = [load(filemill_0_40(:)) ; load(filemill_1_40(:)) ; load(filemill_2_40(:)) ; load(filemill_3_40(:))];

Emill40 = datamill2(:,1);
Cvmill40 = datamill2(:,2);
absMmill40 = datamill2(:,5);
absXmill40 = datamill2(:,6);

% N=60, #MCC = 1e6
filemill_0_60 = 'millnum_60_0.txt';
filemill_1_60 = 'millnum_60_1.txt';
filemill_2_60 = 'millnum_60_2.txt';
filemill_3_60 = 'millnum_60_3.txt';

datamill3 = [load(filemill_0_60(:)) ; load(filemill_1_60(:)) ; load(filemill_2_60(:)) ; load(filemill_3_60(:))];

Emill60 = datamill3(:,1);
Cvmill60 = datamill3(:,2);
absMmill60 = datamill3(:,5);
absXmill60 = datamill3(:,6);

% N=80, #MCC = 1e6
filemill_0_80 = 'millnum_80_0.txt';
filemill_1_80 = 'millnum_80_1.txt';
filemill_2_80 = 'millnum_80_2.txt';
filemill_3_80 = 'millnum_80_3.txt';

datamill4 = [load(filemill_0_80(:)) ; load(filemill_1_80(:)) ; load(filemill_2_80(:)) ; load(filemill_3_80(:))];

Emill80 = datamill4(:,1);
Cvmill80 = datamill4(:,2);
absMmill80 = datamill4(:,5);
absXmill80 = datamill4(:,6);


% figure(2)
% plot(T,Cvmill20,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<Cv/k>','fontsize', 18)
% title('<Cv/k> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% hold('on')
% plot(T,Cvmill40,'co-')
% plot(T,Cvmill60,'go-')
% plot(T,Cvmill80,'mo-')
% legend('N = 20', 'N = 40', 'N = 60', 'N = 80','location','northwest')
% 
% figure(3)
% plot(T,Emill20,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<E/J>','fontsize', 18)
% title('<E/J> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% hold('on')
% plot(T,Emill40,'co-')
% plot(T,Emill60,'go-')
% plot(T,Emill80,'mo-')
% legend('N = 20', 'N = 40', 'N = 60', 'N = 80','location','northwest')
% 
% figure(4)
% plot(T,absMmill20,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<|M|>','fontsize', 18)
% title('<|M|> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% hold('on')
% plot(T,absMmill40,'co-')
% plot(T,absMmill60,'go-')
% plot(T,absMmill80,'mo-')
% legend('N = 20', 'N = 40', 'N = 60', 'N = 80','location','northeast')
% 
% figure(5)
% plot(T,absXmill20,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<|X|J>','fontsize', 18)
% title('<|X|J> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% hold('on')
% plot(T,absXmill40,'co-')
% plot(T,absXmill60,'go-')
% plot(T,absXmill80,'mo-')
% legend('N = 20', 'N = 40', 'N = 60', 'N = 80','location','northwest')

% N=100, #MCC = 1e7
filetmill_0_100 = 'tmillnum_100_0.txt';
filetmill_1_100 = 'tmillnum_100_1.txt';
filetmill_2_100 = 'tmillnum_100_2.txt';
filetmill_3_100 = 'tmillnum_100_3.txt';

datatmill = [load(filetmill_0_100(:)) ; load(filetmill_1_100(:)) ; load(filetmill_2_100(:)) ; load(filetmill_3_100(:))];

E100 = datatmill(:,1);
Cv100 = datatmill(:,2);
absM100 = datatmill(:,5);
absX100 = datatmill(:,6);

% figure(6)
% plot(T,Cv100,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<Cv/k>','fontsize', 18)
% title('<Cv/k> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% 
% figure(7)
% plot(T,E100,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<E/J>','fontsize', 18)
% title('<E/J> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% 
% figure(8)
% plot(T,absM100,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<|M|>','fontsize', 18)
% title('<|M|> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% 
% figure(9)
% plot(T,absX100,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<|X|J>','fontsize', 18)
% title('<|X|J> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% 
% figure(10)
% plot(T,Cvmill20,'ko-')
% xlabel('kT/J','fontsize', 18)
% ylabel('<Cv/k>','fontsize', 18)
% title('<Cv/k> around critical temperature','fontsize', 18)
% set(gca,'FontSize',15)
% hold('on')
% plot(T,Cvmill40,'co-')
% plot(T,Cvmill60,'go-')
% plot(T,Cvmill80,'mo-')
% plot(T,Cv100,'bo-')
% legend('N = 20', 'N = 40', 'N = 60', 'N = 80','N=100','location','northwest')

figure(11)
plot(T,absMmill20,'ko-')
xlabel('kT/J','fontsize', 18)
ylabel('<|M|>','fontsize', 18)
title('<|M|> around critical temperature','fontsize', 18)
set(gca,'FontSize',15)
hold('on')
plot(T,absMmill40,'co-')
plot(T,absMmill60,'go-')
plot(T,absMmill80,'mo-')
plot(T,absM100,'bo-')
legend('N = 20', 'N = 40', 'N = 60', 'N = 80','N = 100','location','northeast')


