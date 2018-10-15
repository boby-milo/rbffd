function [u,err,tim,x,dx,n,N,W] = BSeuCall2Dbasket_RBFFDadap_phs(Nx,g,p,d,nm,M,Kmul)
%% 2D EU Call RBF-FD with BDF2

funname = 'BSeuCall2Dbasket_RBFFDadap_phs';
datafilename = [funname,'___Nx',num2str(Nx),'_g',num2str(g),'_p',num2str(p),'_d',num2str(d),'_nm',num2str(nm),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

names = dir(['./Data/',datafilename]);

if ~isempty(names)
    cd('./Data')
    load(datafilename);
    cd('..')
else
    load('UrefEU.mat')

    tic
    %% Model
    r=0.03;
    sig1=0.15;
    sig2=0.15;
    rho=0.5; %rho=[1 0.5; 0.5 1];

    T=0.2;
    K=1;

    %% Grid
    Kx=1/Kmul;

    % Nx=100;
    i=1:Nx;
    Ki=2*Kx;
    S=1;

    % g=5; %tune this! 1,2,3,4,5

    c=2*Ki/g

    dxi=(1/Nx)*(asinh((S-Ki)/c)-asinh(-Ki/c));
    xi=asinh(-Ki/c)+i*dxi;
    x=[0, Ki+c*sinh(xi)];
    y=zeros(numel(x),1);
    Kind=-x(x<=Ki)+Ki;

    xvec=[]; yvec=[];
    for ii=1:numel(x)
        xl=linspace(0,x(ii),ii);
        yl=linspace(x(ii),0,ii);

        xvec=[xvec,xl];
        yvec=[yvec,yl];
    end

    N=numel(xvec);
    ind=1:N;

    indcf=1;
    indff=(N-numel(x)+1):N;
    indin=ind; indin([indff,indcf])=[];

    dx = 1/sqrt(N);

    % M=100000;
    dt=T/(M-1);
    t=T:-dt:0;
    
    s = [xvec' yvec'];
    
    indc1 = findKNearestNeighbors(s,s(973,:),75);
    indc2 = findKNearestNeighbors(s,s(387,:),75);
    indc3 = findKNearestNeighbors(s,s(253,:),75);
    indc4 = findKNearestNeighbors(s,s(948,:),75);
    
    xvec = 8*xvec; yvec = 8*yvec;
    figure(2)
    clf
    plot(xvec,yvec,'k.')
    hold on
    figure(2); hold on; plot(xvec(indc1(1)),yvec(indc1(1)),'bx'); plot(xvec(indc1),yvec(indc1),'o'); axis equal; axis tight
    figure(2); hold on; plot(xvec(indc2(1)),yvec(indc2(1)),'bx'); plot(xvec(indc2),yvec(indc2),'o'); axis equal; axis tight
    figure(2); hold on; plot(xvec(indc3(1)),yvec(indc3(1)),'bx'); plot(xvec(indc3),yvec(indc3),'o'); axis equal; axis tight
    figure(2); hold on; plot(xvec(indc4(1)),yvec(indc4(1)),'bx'); plot(xvec(indc4),yvec(indc4),'o'); axis equal; axis tight

    %% Initial condition
    u0=max((1/2)*(xvec+yvec)-Kx,zeros(1,length(xvec)));
    u=u0';

    % figure(1)
    % clf
    % plot(xvec,yvec,'.')
    % hold on
    % plot(xvec(indff),yvec(indff),'*')
    % plot(xvec(indcf),yvec(indcf),'^')
    % plot(xvec(indin),yvec(indin),'o')
    % axis equal
    % axis tight
    % hold off

    % figure(2)
    % tri = delaunay(xvec',yvec');
    % trisurf(tri, xvec', yvec', u);
    % shading interp
    % colorbar
    % view(2)
    % axis vis3d
    % % axis equal
    % axis tight

    %% RBF
    phi='phs';

    dim = 2; %problem dimension

    m = nchoosek(p+dim, p); %number of polynomial terms;
    n = round(nm*m);
    s = [xvec' yvec'];

    parallel = 0;
    [W,~] = b(r,sig1,sig2,rho,s,x,y,N,n,m,p,indin,phi,d,'reg',parallel,Kx);

    %% Integration
    I = speye(size(W));
    % BDF-1
    u1 = u;
    A = I-dt*W;

    b = u1;
    b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);

    u = A\b;
    u = max(u,zeros(size(u)));

    % BDF-2
    A = I-(2/3)*dt*W;
    rcm = symrcm(A);
    A = A(rcm,rcm);
    [L1, U1] = lu(A);
    for ii = 3:M
        u2 = u1;
        u1 = u;
        b = (4/3)*u1-(1/3)*u2;
        b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);

        u(rcm) = L1\b(rcm);
        u(rcm) = U1\u(rcm);

        u = max(u,zeros(size(u)));
    end
    tim = toc;


    %% Error
    indreg = [];
    for ii = 1:length(xvec)
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if xvec(ii)>=2/3*Kx && xvec(ii)<=4/3*Kx && yvec(ii)>=2/3*Kx && yvec(ii)<=4/3*Kx
            indreg = [indreg ii];
        end
    end
    xvec=xvec(indreg);
    yvec=yvec(indreg);
    u=u(indreg);

    xvec = K*Kmul*xvec;
    yvec = K*Kmul*yvec;
    dx = K*Kmul*dx;
    u = K*Kmul*u;


    x = [xvec' yvec'];


    uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

    err = uinterp'-u;

    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'n', 'N', 'W', 'Nx', 'g', 'p', 'd', 'nm', 'M', 'Kmul')
    cd('..')
end
end
