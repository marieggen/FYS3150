clear all
close all 
clc

file = 'P_E_T2p4.txt';
P = load(file(:));
E = linspace(-800,800,1600);
plot(E,P,'r')


