ConvergenceTest_phs
N = [40, 80, 120, 160, 200, 240, 280, 320];
M = 4*N;
for ii = 1:hN
    ii
    Nin = N(ii);
    Min = M(ii);
    for jj = 1:numel(p)
        for kk = 1:numel(d)
            for ll = 1:numel(nm)
                for mm = 1:numel(g)
                    [u,err,tim,x,dx,Ntot,W] = feval(leg{4}, N(ii), g(mm), p(jj), d(kk), nm(ll), M(ii), Kmul);
                    
                    cd ('./Data')
                    save([basename,leg{4},'___N',num2str(Nin),'_p',num2str(p(jj)),'_d',num2str(d(kk)),'_nm',num2str(nm(ll)),'_g',num2str(g(mm)),'.mat'], 'u', 'err', 'tim', 'x', 'dx', 'Ntot','W', 'Nin', 'Min', 'Kmul')
                    cd('..')
                end
            end
        end
    end
end
