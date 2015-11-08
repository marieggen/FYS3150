clear all
close all
clc

files = ['ana_b.txt';'num_b.txt'];
files_test = ['ana_test.txt';'num_test.txt'];
files_mcc = ['ana_mcc.txt';'num_mcc.txt'];

ana = load(files(1,:));
num = load(files(2,:));

anaE = ana(:,1);
anaCv = ana(:,2);
anaM = ana(:,3);
anaX = ana(:,4);
anaAbsM = ana(:,5);
anaAbsX = ana(:,6);
T = ana(:,7);
% MCC = ana(:,8)

numE = num(:,1);
numCv = num(:,2);
numM = num(:,3);
numX = num(:,4);
numAbsM = num(:,5);
numAbsX = num(:,6);

%Find relative error in data
eE = (numE-anaE)/anaE
eCv = (numCv-anaCv)/anaCv
eAbsM = (numAbsM-anaAbsM)/anaAbsM
eAbsX = (numAbsX-anaAbsX)/anaAbsX


% figure(1)
% plot(T,anaE,'b-')
% hold('on')
% plot(T,numE,'r-')
% title('Mean energy')
% xlabel('Tk/J')
% ylabel('<E/J>')
% legend('Analytical','Numerical','Location','northwest')
% 
% figure(2)
% plot(T,anaCv,'b-')
% hold('on')
% plot(T,numCv,'r-')
% title('Mean specific heat capacity')
% xlabel('Tk/J')
% ylabel('<Cv/k>')
% legend('Analytical','Numerical','Location','northwest')
% 
% figure(3)
% plot(T,anaM,'b-')
% hold('on')
% plot(T,numM,'r-')
% title('Mean magnetization')
% xlabel('Tk/J')
% ylabel('<M>')
% legend('Analytical','Numerical')
% 
% figure(4)
% plot(T,anaAbsM,'b-')
% hold('on')
% plot(T,numAbsM,'r-')
% title('Mean absolute magnetization')
% xlabel('Tk/J')
% ylabel('<|M|>')
% legend('Analytical','Numerical')
% 
% figure(5)
% plot(T,anaAbsX,'b-')
% hold('on')
% plot(T,numAbsX,'r-')
% title('Mean absolute susceptibility')
% xlabel('Tk/J')
% ylabel('<|X|>')
% legend('Analytical','Numerical','Location','northwest')












