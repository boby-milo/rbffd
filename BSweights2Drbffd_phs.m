function W = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,N,n,m,p,indin,phi,ep,adap,parallel)
% Constructs a BS differentiation matrix W from 2D grid s,
%stencil size n and indices indin.

if nargin == 12
    adap = 'reg';
    parallel = 0;
elseif nargin == 11
    parallel = 0;
end

% if parallel
%     p = gcp();
%     argfor = p.NumWorkers;
% else
%     p = gcp('nocreate');
% %     delete(p);
%     argfor = 0;
% end

% Weights
indc = findKNearestNeighbors(s,s,n);

iind = repmat(indin,n,1); iind = iind(:); %n*N
jind = transpose(indc(indin,:)); jind = jind(:);%n*N
Wval = zeros(n,numel(indin));  %n*N

% parfor (ii = indin, argfor)
for ii = indin
    sc = s(ii,:); xc = sc(:,1); yc = sc(:,2);
    se = s(indc(ii,:),:);
    
    Rc = xcdist(se,se,1);
    
    if adap == 'min'
        H = Rc(:,:,1);
        hmin = min(H(H>0));
        ep = gamma/hmin;
        
    elseif adap == 'avg'
        H = Rc(:,:,1);
        havg = mean(H(H>0));
        ep = gamma/havg;
    end
    
    A = RBFmat(phi,ep,Rc,'0',1);
    
    Ax = RBFmat(phi,ep,Rc,'1',1);
    Ay = RBFmat(phi,ep,Rc,'1',2);
    
    Axx = RBFmat(phi,ep,Rc,'2',1);
    Ayy = RBFmat(phi,ep,Rc,'2',2);
    Axy = RBFmat(phi,ep,Rc,'m2',1:2);
    
    [P,lP] = vander2D(sc,se,n,p,r,sig1,sig2,rho);
    
    
    Ac = [A, P;
          P', zeros(m,m)];
    
    lA = transpose(-r*A(1,:)...
        +r*xc'.*Ax(1,:)+r*yc'.*Ay(1,:)...
        +0.5*sig1^2*xc'.^2.*Axx(1,:)...
        +0.5*sig2^2*yc'.^2.*Ayy(1,:)...
        +rho*sig1*sig2*xc'.*yc'.*Axy(1,:));

    
    
   lc = [lA; lP];
    
    wc = Ac\lc;
    Wval(:,ii-1) = wc(1:n);
end

Wval = Wval(:);
W = sparse(iind,jind,Wval,N,N);
end

function [P, lP] = vander2D(sc,s,n,p,r,sig1,sig2,rho)
x = s(:,1);
y = s(:,2);

xc = sc(1);
yc = sc(2);

dim = 2;
m = nchoosek(p+dim, p); %number of polynomial terms;

P = zeros(n,m);
lP = zeros(m,1);

kk = 0;
for ii = 0:p
    for jj = 0:ii
        kk = kk+1;
        
        i = ii-jj;
        j = jj;
        
        P(:,kk) = x.^(i) .* y.^(j);
        
        lP(kk) = xc^(i)*yc^(j) * (-r +r*(i+j) + 0.5*(i-1)*i*sig1^2 + 0.5*(j-1)*j*sig2^2 +i*j*rho*sig1*sig2);
        
    end
end
end