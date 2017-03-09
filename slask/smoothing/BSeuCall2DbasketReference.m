clc
clear all
close all

Kmul = 4;

Nx = 760 
Nx1 = round(0.9*760);

M = 10^5;

[u1,~,tim1,x1,dx1,n1,N1,W1] = BSeuCall2Dbasket_FD(Nx1,M,Kmul)
[u,~,tim,x,dx,n,N,W] = BSeuCall2Dbasket_FD(Nx,M,Kmul)