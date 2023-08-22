clear all
close all
clc

% L2 = [0.9876 0.2492 0.0625 0.0156 0.0039];
% H1 = [3.4157 1.726 0.8653 0.4329 0.2165];
% h = [1 0.5 0.25 0.125 0.0625];

L2 = [1.0806 0.2729 0.0684 0.0171 0.0043];
H1 = [3.4254 1.7272 0.8654 0.4329 0.2165];
h = [1 0.5 0.25 0.125 0.0625];
dof = [4 7 13 25 49];

qL2 = [0.0333 0.0042 5.2083e-04];
qH1 = [0.2582 0.0645 0.0161];
qh = [1 0.5 0.25];
qdof = [4 7 13];

% loglog(h,L2,'*--b')
% hold on
% loglog(h,H1,'*--r')
% polyfit(log(h),log(L2),1)
% polyfit(log(h),log(H1),1)
% 
% loglog(qh,qL2,'*-b')
% hold on
% loglog(qh,qH1,'*-r')
% polyfit(log(qh),log(qL2),1)
% polyfit(log(qh),log(qH1),1)

loglog(dof,L2,'*--b')
hold on
loglog(dof,H1,'*--r')
polyfit(log(dof),log(L2),1)
polyfit(log(dof),log(H1),1)

loglog(qdof,qL2,'*-b')
hold on
loglog(qdof,qH1,'*-r')
polyfit(log(qdof),log(qL2),1)
polyfit(log(qdof),log(qH1),1)
