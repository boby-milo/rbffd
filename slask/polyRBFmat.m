function A = polyRBFmat(phi,ep,s,nprime,dimprime,p)
% s - evaluation points [N x dim]
% ep - shape parameter
% phi - radial basis function [string]
% nprime - order of derivative [string]
% p - order of polynomials [integer]

N = size(s,1);
dim = size(s,2);

R = xcdist(s,s,1);
A = feval(phi,ep,R,nprime,dimprime);

m = nchoosek(p+dim, p); %number of polynomial terms;
P = zeros(N,m);

switch dim
    case 1
        switch nprime
            case '0'
                ll = 0;
                for kk = 0:p
                    ll = ll+1;
                    P(:,ll) = s(:,1).^(kk);
                end
                
            case '1'
                ll = 0;
                for kk = 0:p
                    ll = ll+1;
                    P(:,ll) = kk*s(:,1).^(kk-1);
                end
                
                
            case '2'
                ll = 0;
                for kk = 0:p
                    ll = ll+1;
                    P(:,ll) = (kk-1)*(kk)*s(:,1).^(kk-2);
                end
        end
        
        
        
    case 2
        switch nprime
            case '0'
                ll = 0;
                for jj = 0:p
                    for kk = 0:p
                        if jj+kk <= p
                            ll = ll+1;
                            P(:,ll) = s(:,1).^(kk) .* s(:,2).^(jj);
                        end
                    end
                end
                
            case '1'
                if dimprime == 1
                    ll = 0;
                    for jj = 0:p
                        for kk = 0:p
                            if jj+kk <= p
                                ll = ll+1;
                                P(:,ll) = kk*s(:,1).^(kk-1) .* s(:,2).^(jj);
                            end
                        end
                    end
                    
                elseif dimprime == 2
                    ll = 0;
                    for jj = 0:p
                        for kk = 0:p
                            if jj+kk <= p
                                ll = ll+1;
                                P(:,ll) = s(:,1).^(kk) .* (jj*s(:,2).^(jj-1));
                            end
                        end
                    end
                end
                
            case '2'
                if dimprime == 1
                    ll = 0;
                    for jj = 0:p
                        for kk = 0:p
                            if jj+kk <= p
                                ll = ll+1;
                                P(:,ll) = (kk-1)*kk*s(:,1).^(kk-2) .* s(:,2).^(jj);
                            end
                        end
                    end
                    
                elseif dimprime == 2
                    ll = 0;
                    for jj = 0:p
                        for kk = 0:p
                            if jj+kk <= p
                                ll = ll+1;
                                P(:,ll) = s(:,1).^(kk) .* ((jj-1)*jj*s(:,2).^(jj-2));
                            end
                        end
                    end
                end
                
            case 'm2'
                ll = 0;
                for jj = 0:p
                    for kk = 0:p
                        if jj+kk <= p
                            ll = ll+1;
                            P(:,ll) = kk*s(:,1).^(kk-1) .* (jj*s(:,2).^(jj-1));
                        end
                    end
                end
        end
        
        
        
    case 3
        switch nprime
            case '0'
                ll = 0;
                for ii = 0:p
                    for jj = 0:p
                        for kk = 0:p
                            if ii+jj+kk <= p
                                ll = ll+1;
                                P(:,ll) = s(:,1).^(kk) .* s(:,2).^(jj) .* s(:,3).^(ii);
                            end
                        end
                    end
                end
                
            case '1'
                if dimprime == 1
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = kk*s(:,1).^(kk-1) .* s(:,2).^(jj) .* s(:,3).^(ii);
                                end
                            end
                        end
                    end
                    
                elseif dimprime == 2
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = s(:,1).^(kk) .* (jj*s(:,2).^(jj-1)) .* s(:,3).^(ii);
                                end
                            end
                        end
                    end
                    
                elseif dimprime == 3
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = s(:,1).^(kk) .* s(:,2).^(jj) .* (ii*s(:,3).^(ii-1));
                                end
                            end
                        end
                    end
                end
                
                
            case '2'
                if dimprime == 1
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = (kk-1)*kk*s(:,1).^(kk-2) .* s(:,2).^(jj) .* s(:,3).^(ii);
                                end
                            end
                        end
                    end
                    
                elseif dimprime == 2
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = s(:,1).^(kk) .* ((jj-1)*jj*s(:,2).^(jj-2)) .* s(:,3).^(ii);
                                end
                            end
                        end
                    end
                    
                elseif dimprime == 3
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = s(:,1).^(kk) .* s(:,2).^(jj) .* ((ii-1)*ii*s(:,3).^(ii-2));
                                end
                            end
                        end
                    end
                end
                
            case 'm2'
                if norm(dimprime - [1,2]) == 0
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = kk*s(:,1).^(kk-1) .* (jj*s(:,2).^(jj-1)) .* s(:,3).^(ii);
                                end
                            end
                        end
                    end
                    
                elseif norm(dimprime - [1,3]) == 0
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = kk*s(:,1).^(kk-1) .* s(:,2).^(jj) .* (ii*s(:,3).^(ii-1));
                                end
                            end
                        end
                    end
                    
                elseif norm(dimprime - [2,3]) == 0
                    ll = 0;
                    for ii = 0:p
                        for jj = 0:p
                            for kk = 0:p
                                if ii+jj+kk <= p
                                    ll = ll+1;
                                    P(:,ll) = s(:,1).^(kk) .* (jj*s(:,2).^(jj-1)) .* (ii*s(:,3).^(ii-1));
                                end
                            end
                        end
                    end
                    
                end
        end
        
        A = [A, P;
            P', zeros(m,m)];
        
end