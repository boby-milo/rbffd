function [ep] = epsopt(h,n,D)
%Returns optimal shape parameter for regular stencils

if D==1
    if n==3
        p1 = 0.0005057;
        p2 = 0.3446;
        ep = 4*p1*(1./h) + p1+p2;
    elseif n==5
        
    elseif n==7
        p1 = 0.01314;
        p2 = 7.589;
        ep = 4*p1*(1./h) + p1+p2;
    elseif n==9
        
    end
end

