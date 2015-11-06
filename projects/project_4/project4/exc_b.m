clear all
close all

files = ['ana.txt';'num.txt'];%Temp = 1.0;
ana = load(files(1,:));
num = load(files(2,:));

anaE = ana(:,1);
anaCv = ana(:,2);
anaM = ana(:,3);
anaX = ana(:,4);
anaAbsM = ana(:,5);
anaAbsX = ana(:,6);
T = ana(:,7);

numE = num(:,1);
numCv = num(:,2);
numM = num(:,3);
numX = num(:,4);
numAbsM = num(:,5);
numAbsX = num(:,6);

figure(1)
plot(T,anaE)
hold('on')
plot(T,numE)

figure(2)
plot(T,anaCv)
hold('on')
plot(T,numCv)

figure(3)
plot(T,anaM)
hold('on')
plot(T,numM)

figure(4)
plot(T,anaX)
hold('on')
plot(T,numX)

figure(5)
plot(T,anaAbsM)
hold('on')
plot(T,numAbsM)

figure(6)
plot(T,anaAbsX)
hold('on')
plot(T,numAbsX)












