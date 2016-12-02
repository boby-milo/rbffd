function W = BSweights1Drbffd3poly(r,sig,x,N,n,indin,phi,ep,parallel)
% Constructs a BS differentiation matrix W from 1D grid x,
%stencil size n.

if parallel
    p = gcp();
    argfor = p.NumWorkers;
else
    p = gcp('nocreate');
    delete(p);
    argfor = 0;
end

m=round((n-1)/2);

% Weights
iind=repmat(indin,n,1); iind=iind(:);
jind=zeros(n,N-2);
Wval=zeros(n,N-2);

for ii=2:m
    xc=x(ii);
    indc=1:n;
    a=ii;
    
    wc = RBFelements(x,xc,n,indc,phi,ep,r,sig,a);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

parfor (ii=(m+1):(N-m), argfor)
    xc=x(ii);
    indc=ii-m:ii+m;
    a=m+1;
    
    wc = RBFelements(x,xc,n,indc,phi,ep,r,sig,a);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

for ii=(N-m+1):(N-1)
    xc=x(ii);
    indc=N-n+1:N;
    a=ii-N+n;
    
    wc = RBFelements(x,xc,n,indc,phi,ep,r,sig,a);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

jind=jind(:);
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);
end

function wc = RBFelements(x,xc,n,indc,phi,ep,r,sig,a)
lc=zeros(n+1,1);

Rc=xcdist(x(indc),x(indc),1);

%     H=Rc(:,:,1);
%     hmin=min(min(H(H>0)));
%     ep=gamma/hmin;

A=RBFmat(phi,ep,Rc,'0',1);
Ax=RBFmat(phi,ep,Rc,'1',1);
Axx=RBFmat(phi,ep,Rc,'2',1);

o = ones(1,n);
o1 = transpose(x(indc));
o2 = transpose(x(indc).^2);
Ac=[A, transpose(o), transpose(o1), transpose(o2);
    o, zeros(1,3);
    o1, zeros(1,3);
    o2, zeros(1,3);];

lc(1:n,1) = transpose(-r*A(a,:)+r*xc.*Ax(a,:)+0.5*xc.^2.*sig^2.*Axx(a,:));
lc(n+1,1) = -r;
lc(n+2,1) = 0;
lc(n+3,1) = r*xc.^2 + sig^2;

wc=Ac\lc;
wc=wc(1:end-3);
end