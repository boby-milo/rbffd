function W = BSweights2Drbffd(r,sig1,sig2,rho,s,N,n,indin,phi,ep,adap)
% Constructs a BS differentiation matrix W from 2D grid s,
%stencil size n and indices indin.

% Weights
indc=findKNearestNeighbors(s,s,n);

iind=repmat(indin,n,1); iind=iind(:); %n*N
jind=transpose(indc(indin,:)); jind=jind(:);%n*N
Wval=zeros(n,numel(indin));  %n*N

lc=zeros(n+1,1);

for ii=indin
    sc=s(ii,:); xc=sc(:,1); yc=sc(:,2);
    se=s(indc(ii,:),:);
    
    Rc=xcdist(se,se,1);
    
    if adap == 'min'
        H=Rc(:,:,1);
        hmin=min(H(H>0));
        ep=gamma/hmin;
        
    elseif adap == 'avg'
        H=Rc(:,:,1);
        havg=mean(H(H>0));
        ep=gamma/havg;
        
    end
    
    A=RBFmat(phi,ep,Rc,'0',1);
    
    Ax=RBFmat(phi,ep,Rc,'1',1);
    Ay=RBFmat(phi,ep,Rc,'1',2);
    
    Axx=RBFmat(phi,ep,Rc,'2',1);
    Ayy=RBFmat(phi,ep,Rc,'2',2);
    Axy=RBFmat(phi,ep,Rc,'m2',1:2);
    
    o = ones(1,n);
    Ac = [A, transpose(o);
        o, 0];
    
    lc(1:n,1) = transpose(-r*A(1,:)...
        +r*xc'.*Ax(1,:)+r*yc'.*Ay(1,:)...
        +0.5*sig1^2*xc'.^2.*Axx(1,:)...
        +0.5*sig2^2*yc'.^2.*Ayy(1,:)...
        +rho*sig1*sig2*xc'.*yc'.*Axy(1,:));
    
    lc(n+1,1) = -r;
    
    wc = Ac\lc;
    Wval(:,ii-1) = wc(1:end-1);
end

Wval = Wval(:);
W = sparse(iind,jind,Wval,N,N);
end
