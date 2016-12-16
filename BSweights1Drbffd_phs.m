function W = BSweights1Drbffd_phs(r,sig,x,N,n,m,indin,phi,ep,parallel)
% Constructs a BS differentiation matrix W from 1D grid x,
%stencil size n.

if nargin == 8
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

l=round((n-1)/2); %stencil distance;

% Weights
iind=repmat(indin,n,1); iind=iind(:);
jind=zeros(n,N-2);
Wval=zeros(n,N-2);

for ii=2:l
    xc=x(ii);
    indc=1:n;
    a=ii;
    
    wc = RBFelements(x,xc,n,m,indc,phi,ep,r,sig,a);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

% parfor (ii=(l+1):(N-l), argfor)
for ii=(l+1):(N-l)
    xc=x(ii);
    indc=ii-l:ii+l;
    a=l+1;
    
    wc = RBFelements(x,xc,n,m,indc,phi,ep,r,sig,a);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

for ii=(N-l+1):(N-1)
    xc=x(ii);
    indc=N-n+1:N;
    a=ii-N+n;
    
    wc = RBFelements(x,xc,n,m,indc,phi,ep,r,sig,a);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

jind=jind(:);
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);
end

function wc = RBFelements(x,xc,n,m,indc,phi,ep,r,sig,a)

Rc=xcdist(x(indc),x(indc),1);

%     H=Rc(:,:,1);
%     hmin=min(min(H(H>0)));
%     ep=gamma/hmin;

A=RBFmat(phi,ep,Rc,'0',1);
Ax=RBFmat(phi,ep,Rc,'1',1);
Axx=RBFmat(phi,ep,Rc,'2',1);

P = vander(x(indc));
P = fliplr(P);
P = P (:,1:m);


Ac=[A, P;
    P', zeros(m,m)];


lphi = transpose(-r*A(a,:)+r*xc.*Ax(a,:)+0.5*xc.^2.*sig^2.*Axx(a,:));

i = transpose(1:m);
lp=((i-2)*r + 0.5*(i-2).*(i-1)*sig^2).*xc.^(i-1);


lc = [lphi;lp];


wc=Ac\lc;
wc=wc(1:n);
end