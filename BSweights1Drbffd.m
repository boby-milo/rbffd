function W = BSweights1Drbffd(r,sig,x,N,n,indin,phi,ep)
% Constructs a BS differentiation matrix W from 1D grid x,
%stencil size n.

m=round((n-1)/2);

% Weights
iind=repmat(indin,n,1); iind=iind(:);
jind=zeros(n,N-2);
Wval=zeros(n,N-2);

for ii=2:m
    xc=x(ii);
    indc=1:n;
    a=ii;
    
    [Wval, jind] = RBFelements(ii,x,xc,n,indc,phi,ep,r,sig,jind,Wval,a);
end

for ii=(m+1):(N-m)
    xc=x(ii);
    indc=ii-m:ii+m;
    a=m+1;
    
    [Wval, jind] = RBFelements(ii,x,xc,n,indc,phi,ep,r,sig,jind,Wval,a);
end

for ii=(N-m+1):(N-1)
    xc=x(ii);
    indc=N-n+1:N;
    a=ii-N+n;
    
    [Wval, jind] = RBFelements(ii,x,xc,n,indc,phi,ep,r,sig,jind,Wval,a);
end

jind=jind(:);
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);
end

function [Wval, jind] = RBFelements(ii,x,xc,n,indc,phi,ep,r,sig,jind,Wval,a)
lc=zeros(n+1,1);

Rc=xcdist(x(indc),x(indc),1);

%     H=Rc(:,:,1);
%     hmin=min(min(H(H>0)));
%     ep=gamma/hmin;

A=RBFmat(phi,ep,Rc,'0',1);
Ax=RBFmat(phi,ep,Rc,'1',1);
Axx=RBFmat(phi,ep,Rc,'2',1);

o=ones(1,n);
Ac=[A, transpose(o);
    o, 0];

lc(1:n,1)=transpose(-r*A(a,:)+r*xc.*Ax(a,:)+0.5*xc.^2.*sig^2.*Axx(a,:));
lc(n+1,1)=-r;

wc=Ac\lc;

Wval(:,ii-1)=wc(1:end-1);

jind(:,ii-1)=indc';

end
