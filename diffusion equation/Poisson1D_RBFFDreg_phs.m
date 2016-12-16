function [u,err,tim,x,dx,N,W] = Poisson1D_RBFFDreg_phs(N,p)
%% 1D European Call RBF-FD
% 2016-02-06

tic
%% Grid
x = transpose(linspace(-1, 1, N));
dx = x(2) - x(1);
indin = 2 : N-1;

%% Parameters
uexact = sin(10*x); 

f = -100*sin(10*x); %force term

u0 = sin(10*x(1)) %boundary
uend = sin(10*x(end)) %boundary

b = f; b(1) = u0; b(end) = uend;

%% RBF
phi = 'phs';
d = 3;

% p = 1;
m = p + 1; %number of polynomial terms;

n = 2 * m;
if mod(n, 2)
    n
else
    n = n + 1
end
%stencil size

l = round((n-1)/2);

iind=repmat(indin,n,1); iind=iind(:);
jind=zeros(n,N-2);
Wval=zeros(n,N-2);

for ii=2:l
    xc=x(ii);
    indc=1:n;
    jjc=ii;
    
%     figure(1)
%     plot(x,zeros(size(x)),'k.');
%     hold on
%     plot(xc,zeros(size(xc)),'rx');
%     plot(x(indc),zeros(size(x(indc))),'go');
%     drawnow;
%     pause(0.01)
%     clf
    
    wc = RBFelements(x,xc,n,m,indc,phi,d,jjc);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

% parfor (ii=(l+1):(N-l), argfor)
for ii=(l+1):(N-l)
    xc=x(ii);
    indc=ii-l:ii+l;
    jjc=l+1;
    
%     figure(1)
%     plot(x,zeros(size(x)),'k.');
%     hold on
%     plot(xc,zeros(size(xc)),'rx');
%     plot(x(indc),zeros(size(x(indc))),'go');
%     drawnow;
%     pause(0.01)
%     clf
    
    wc = RBFelements(x,xc,n,m,indc,phi,d,jjc);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

for ii=(N-l+1):(N-1)
    xc=x(ii);
    indc=N-n+1:N;
    jjc=ii-N+n;
    
%     figure(1)
%     plot(x,zeros(size(x)),'k.');
%     hold on
%     plot(xc,zeros(size(xc)),'rx');
%     plot(x(indc),zeros(size(x(indc))),'go');
%     drawnow;
%     pause(0.01)
%     clf
    
    wc = RBFelements(x,xc,n,m,indc,phi,d,jjc);
    Wval(:,ii-1)=wc;
    jind(:,ii-1)=indc';
end

jind=jind(:);
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);

W(1,1) = 1;
W(end,end) = 1;

tim = toc;

u = W\b;

err = uexact - u;

end


function wc = RBFelements(x,xc,n,m,indc,phi,d,jjc)

Rc=xcdist(x(indc),x(indc),1);

%     H=Rc(:,:,1);
%     hmin=min(min(H(H>0)));
%     d=gamma/hmin;

A = RBFmat(phi,d,Rc,'0',1);
% Ax=RBFmat(phi,d,Rc,'1',1);
Axx = RBFmat(phi,d,Rc,'2',1);

P = vander(x(indc));
P = P (:,1:m);

Ac=[A, P;
    P', zeros(m,m)];

lphi = transpose(Axx(jjc,:));

i = transpose(1:m);
lp=(i-2).*(i-1).*xc.^(i-3);

lc = [lphi;lp];

wc=Ac\lc;
wc=wc(1:n);
end