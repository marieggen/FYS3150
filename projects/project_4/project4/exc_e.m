clear all
close all
clc

file_0 = 'num_test_0.txt';
file_1 = 'num_test_1.txt';
file_2 = 'num_test_2.txt';
file_3 = 'num_test_3.txt';

data1 = load(file_0(:));

data = [load(file_0(:)) ; load(file_1(:)) ; load(file_2(:)) ; load(file_3(:))];

E = data(:,1);
Cv = data(:,2);
absM = data(:,5);
absX = data(:,6);
T = data(:,7);
plot(T,Cv)