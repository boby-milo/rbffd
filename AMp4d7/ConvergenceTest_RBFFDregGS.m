ConvergenceTest_phs

n = 25;

for ii = 1:hN
    
    Nin = N(ii);
    Min = M(ii);
    [u,err,tim,x,dx,ep,Ntot,W] = feval('BSeuCallbasket2D_RBFFDreg', N(ii), n, M(ii), Kmul);
    
end
