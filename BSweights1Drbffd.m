function W = BSweights1Drbffd(r,sig,x,N,n,indin,phi,ep)
% Constructs a BS differentiation matrix W from 1D grid x,
%stencil size n.

m=round((n-1)/2);

Rc=zeros(N,2*n-1,2);

ind=1:n-1;
Rc(ind,1:n+ind-1,:)=xcdist(x(ind),x(1:n+ind-1),1); %#ok<*BDSCI>

ind=n:N-n+1;
Rc(ind,:,:)=xcdist(x(ind),x(ind-n+1:ind+n-1),1);

ind=N-n+2:N;
Rc(ind,ind-N+n:2*n-1,:)=xcdist(x(ind),x(ind-n+1:N),1);

A=RBFmat(phi,ep,Rc,'0',1);
Ax=RBFmat(phi,ep,Rc,'1',1);
Axx=RBFmat(phi,ep,Rc,'2',1);

A = repmat(A(n,:),[N 1]);
A = spdiags(A,-n+1:n-1,N,N);
Ax = repmat(Ax(n,:),[N 1]);
Ax = spdiags(Ax,-n+1:n-1,N,N);
Axx = repmat(Axx(n,:),[N 1]);
Axx = spdiags(Axx,-n+1:n-1,N,N);

% Weights
iind=repmat(indin,n,1); iind=iind(:);
jind=zeros(n,N-2);
Wval=zeros(n,N-2);

lc=zeros(n+1,1);

for ii=2:m
    xc=x(ii);
    indc=1:n;
    o=ones(1,n);
    Ac=[A(indc,indc), transpose(o);
        o, 0];
    lc(1:n,1)=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    lc(n+1,1)=-r;
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    jind(:,ii-1)=indc';
end

for ii=(m+1):(N-m)
    xc=x(ii);
    indc=ii-m:ii+m;
    o=ones(1,n);
    Ac=[A(indc,indc), transpose(o);
        o, 0];
    lc(1:n,1)=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    lc(n+1,1)=-r;
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    jind(:,ii-1)=indc';
end

for ii=(N-m+1):(N-1)
    xc=x(ii);
    indc=N-n+1:N;
    o=ones(1,n);
    Ac=[A(indc,indc), transpose(o);
        o, 0];
    lc(1:n,1)=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    lc(n+1,1)=-r;
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    jind(:,ii-1)=indc';
end

jind=jind(:);
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);
end
