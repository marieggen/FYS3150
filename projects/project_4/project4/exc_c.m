clear all
close all

files_mcc = 'num_mcc.txt';%Temp = 1.0
num = load(files_mcc(:));

numE = num(:,1);
numCv = num(:,2);
numAbsM = num(:,3);
numAbsX = num(:,4);
MCC = num(:,5);

figure(1)
plot(MCC,numE)

% figure(2)
% plot(MCC,numCv)

figure(3)
plot(MCC,numAbsM)

% figure(4)
% plot(MCC,numAbsX)