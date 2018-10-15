ConvergenceTest_phs

n = 25;
% g = 3;
for ii = 1:hN
    
    Nin = N(ii);
    Min = M(ii);
    [u,err,tim,x,dx,ep,Ntot,W] = feval('BSeuCallbasket2D_RBFFDadap', N(ii), g, n, M(ii), Kmul);
    
end
