clear all
close all
clc

files_mcc = 'num_mcc.txt';%Temp = 1.0
num = load(files_mcc(:));

numE = num(:,1)
numCv = num(:,2);
numAbsM = num(:,5);
numAbsX = num(:,6);
MCC = num(:,8);
count = num(:,9);

figure(1)
plot(MCC,numE)

figure(2)
plot(MCC,count)

% figure(2)
% plot(MCC,numCv)

figure(3)
plot(MCC,numAbsM)

% figure(4)
% plot(MCC,numAbsX)