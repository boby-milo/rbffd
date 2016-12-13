clear
clc
close all

Kmul = 4;
N = Kmul*[100, 200, 400, 800]
M = 10000
phi = 'phs'


for ii = 1:numel(N)
    %%
    
    ep = 3
    
    n = 3;
    [u,err,tim,x,dx,N(ii),W] = BSeuCall1D_RBFFDreg_phs(N(ii),n,ep,M);
    errmax = norm(err,Inf)
    
    figure(1)
    loglog(dx,errmax,'r-*')
    hold on
    
    n = 5;
    [u,err,tim,x,dx,N(ii),W] = BSeuCall1D_RBFFDreg_phs(N(ii),n,ep,M);
    errmax = norm(err,Inf)
    
    figure(1)
    loglog(dx,errmax,'g-*')
    hold on
    
    n = 7;
    [u,err,tim,x,dx,N(ii),W] = BSeuCall1D_RBFFDreg_phs(N(ii),n,ep,M);
    errmax = norm(err,Inf)
    figure(1)
    loglog(dx,errmax,'b-*')
    hold on
    
    n = 13;
    [u,err,tim,x,dx,N(ii),W] = BSeuCall1D_RBFFDreg_phs(N(ii),n,ep,M);
    errmax = norm(err,Inf)
      figure(1)
    loglog(dx,errmax,'m-*')
    hold on
    
     n = 15;
    [u,err,tim,x,dx,N(ii),W] = BSeuCall1D_RBFFDreg_phs(N(ii),n,ep,M);
    errmax = norm(err,Inf)
      figure(1)
    loglog(dx,errmax,'y-*')
    hold on
    
%     legend('F1', 'F2', 'F3', 'F4','F5', 'F6', 'F7');
%     title(['M=',num2str(M)]);
%     %     xlim([1e-4, 1e-1]);
%     %     ylim([1e-10, 1e-4]);
    drawnow
end