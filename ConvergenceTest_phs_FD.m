ConvergenceTest_phs

for ii = 1:hN
    ii
    
    Nin = N(ii);
    Min = M(ii);
    [u,err,tim,x,dx,Ntot,W] = feval(leg{1}, N(ii), M(ii), Kmul);
    
    
    cd ('./Data')
    save([basename,leg{1},'___N',num2str(Nin),'.mat'], 'u', 'err', 'tim', 'x', 'dx', 'Ntot','W', 'Nin', 'Min', 'Kmul')
    cd ('..')
end
