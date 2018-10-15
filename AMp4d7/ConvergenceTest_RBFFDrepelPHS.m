ConvergenceTest_phs
warning off
dbstop if error
pp = gcp();
for ii = 1:hN
    Nin = N(ii)
    Min = M(ii)
    for jj = 1:numel(p)
        for kk = 1:numel(d)
            for ll = 1:numel(nm)
                for mm = 1:numel(g)
                    
%                     pp = gcp();
%                     [u,err,tim,x,dx,Ntot,W] = feval('BSeuCall2Dbasket_RBFFDrepel_phs', N(ii), p(jj), d(kk), nm(ll), M(ii), Kmul);
%                     disp('BSeuCall2Dbasket_RBFFDadap_phs')
%                     dx
%                     max(abs(err))
%                     tim
                    
%                     pp = gcp();
%                     [u,err,tim,x,dx,Ntot,W] = feval('BSeuCall2Dbasket_RBFFDrepel_smoothed_phs', N(ii), p(jj), d(kk), nm(ll), M(ii), Kmul);
%                     disp('BSeuCall2Dbasket_RBFFDadap_phs')
%                     dx
%                     max(abs(err))
%                     tim

                    pp = gcp();
                    [u,err,tim,x,dx,Ntot,W] = feval('BSamPut2Dbasket_RBFFDadap_phs', N(ii), g(mm), p(jj), d(kk), nm(ll), M(ii), Kmul);
                    disp('BSamPut2Dbasket_RBFFDadap_phs')
                    dx
                    max(abs(err))
                    tim

                    pp = gcp();
                    [u,err,tim,x,dx,Ntot,W] = feval('BSamPut2Dbasket_RBFFDadap_smoothed_phs', N(ii), g(mm), p(jj), d(kk), nm(ll), M(ii), Kmul);
                    disp('BSamPut2Dbasket_RBFFDadap_smoothed_phs')
                    dx
                    max(abs(err))
                    tim

%                     pp = gcp();
%                     [u,err,tim,x,dx,Ntot,W] = feval('BSamPut2Dbasket_RBFFDadap_gs', N(ii), g(mm), 25, M(ii), Kmul);
%                     disp('BSamPut2Dbasket_RBFFDadap_gs')
%                     dx
%                     max(abs(err))
%                     tim
% 
%                     pp = gcp();
%                     [u,err,tim,x,dx,Ntot,W] = feval('BSamPut2Dbasket_RBFFDreg_gs', N(ii), 25, M(ii), Kmul);
%                     disp('BSamPut2Dbasket_RBFFDreg_gs')
%                     dx
%                     max(abs(err))
%                     tim

                    pp = gcp();
                    [u,err,tim,x,dx,Ntot,W] = feval('BSamPut2Dbasket_RBFFDreg_phs', N(ii), p(jj), d(kk), nm(ll), M(ii), Kmul);
                    disp('BSamPut2Dbasket_RBFFDreg_phs')
                    dx
                    max(abs(err))
                    tim

                    pp = gcp();
                    [u,err,tim,x,dx,Ntot,W] = feval('BSamPut2Dbasket_RBFFDreg_smoothed_phs', N(ii), p(jj), d(kk), nm(ll), M(ii), Kmul);
                    disp('BSamPut2Dbasket_RBFFDreg_smoothed_phs')
                    dx
                    max(abs(err))
                    tim

                end
            end
        end
    end
end
