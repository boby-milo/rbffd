clc
clear all
close all

Kmul = 4;

Nx = 760;
Nx1 = round(0.9*Nx);

M = 10^5;

display(1)
[u1,~,tim1,x1,dx1,n1,N1,W1] = BSeuCall2Dbasket_FD(Nx1,M,Kmul);
display(2)
[u,~,tim,x,dx,n,N,W] = BSeuCall2Dbasket_FD(Nx,M,Kmul);


xvec = x(:,1); yvec = x(:,2);
xvec1 = x1(:,1); yvec1 = x1(:,2);

uint = griddata(xvec,yvec,u,xvec1,yvec1);

err=max(abs(uint-u1))