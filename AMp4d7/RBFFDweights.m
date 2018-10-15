function [W,wcond,hloc] = RBFFDweights(par,s,N,phi,ep,p,n,indin,indc,dim)
% Constructs a SABR differentiation matrix W from 2D grid s,
%stencil size n and indices indin.

ppar = gcp('nocreate');
argfor = ppar.NumWorkers;

m = nchoosek(p+dim, p); %number of polynomial terms

% indc = findKNearestNeighbors(s,s,n);

iind = repmat(indin,n,1); iind = iind(:); %n*N
jind = transpose(indc(indin,:)); jind = jind(:);%n*N
Wval = zeros(n,numel(indin));  %n*N
hloc = zeros(N,1);

parfor (jj = 1:numel(indin), argfor)
    warning off
% for jj = 1:numel(indin)
    ii = indin(jj);
    sc = s(ii,:); %xc = sc(:,1); yc = sc(:,2);
    se = s(indc(ii,:),:);
    
    Rc = xcdist(se,se,1);
    
    H = Rc(:,:,1);
    h(jj) = min(H(H>0));
    
    A = RBFmat(phi,ep,Rc,'0',1);
    Ax = RBFmat(phi,ep,Rc,'1',1);
    Ay = RBFmat(phi,ep,Rc,'1',2);
    Axx = RBFmat(phi,ep,Rc,'2',1);
    Ayy = RBFmat(phi,ep,Rc,'2',2);
    Axy = RBFmat(phi,ep,Rc,'m2',1:2);
    
    
    [P,lP] = vanderRBF(ii,par,sc,se,n,m,p,dim); % defined below!
    
    
    Ac = [A, P;
        P', zeros(m,m)];
    
    lA = transpose(par(ii,1) * A(1,:)...
        + par(ii,2) *     Ax(1,:)...
        + par(ii,3) *     Ay(1,:)...
        + par(ii,4) *     Axx(1,:)...
        + par(ii,5) *     Axy(1,:)...
        + par(ii,6) *     Ayy(1,:));
    
    lc = [lA; lP];
    
    wc = Ac\lc;
    wcond(jj) = cond(Ac);
    
    Wval(:,jj) = wc(1:n);
end

Wval = Wval(:);
W = sparse(iind,jind,Wval,N,N);
hloc(indin) = h;
end

function [P, lP] = vanderRBF(index,par,sc,s,n,m,p,dim)
if dim == 1
    
    
    if (par(index,2) == 0 && par(index,4) == 0)
        x = s(:,2);
        xc = sc(2);
        
        a = par(index,1);
        b = par(index,3);
        c = par(index,6);
        
    elseif (par(index,3) == 0 && par(index,6) == 0)
        x = s(:,1);
        xc = sc(1);
        
        a = par(index,1);
        b = par(index,2);
        c = par(index,4);
    end
    
    %     P = zeros(n,m);
    lP = zeros(m,1);
    
    P = vander(x);
    P = fliplr(P);
    P = P (:,1:m);
    
    for ii = 0:p
        %         if xc == 0
        %             xc1 = 0;
        %             xc2 = 0;
        %         else
        if ii >= 2
            xc1 = xc^(ii-1);
            xc2 = xc^(ii-2);
        elseif ii == 1
            xc1 = xc^(ii-1);
            xc2 = 0;
        else
            xc1 = 0;
            xc2 = 0;
        end
        %         end
        lP(ii+1) = a *      xc^ii...
            + b *       ii*xc1...
            + c * (ii-1)*ii*xc2;
    end
    
    
    
elseif dim == 2
    x = s(:,1);
    y = s(:,2);
    
    xc = sc(1);
    yc = sc(2);
    
    P = zeros(n,m);
    lP = zeros(m,1);
    
    kk = 0;
    for ii = 0:p
        for jj = 0:ii
            kk = kk+1;
            
            i = ii-jj;
            j = jj;
            
            P(:,kk) = x.^(i) .* y.^(j);
            
            if i >= 2
                xc1 = xc^(i-1);
                xc2 = xc^(i-2);
            elseif i == 1
                xc1 = xc^(i-1);
                xc2 = 0;
            else
                xc1 = 0;
                xc2 = 0;
            end
            
            if j >= 2
                yc1 = yc^(j-1);
                yc2 = yc^(j-2);
            elseif j == 1
                yc1 = yc^(j-1);
                yc2 = 0;
            else
                yc1 = 0;
                yc2 = 0;
            end
            
            lP(kk) = par(index,1) *      xc^i *         yc^j...
                + par(index,2) *       i*xc1  *         yc^j...
                + par(index,3) *         xc^i *       j*yc1...
                + par(index,4) * (i-1)*i*xc2  *         yc^j...
                + par(index,5) *       i*xc1  *       j*yc1...
                + par(index,6) *         xc^i * (j-1)*j*yc2;
            
        end
    end
end

end