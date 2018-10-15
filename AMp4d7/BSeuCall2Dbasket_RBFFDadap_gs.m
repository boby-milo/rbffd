function [u,err,tim,x,dx,ep,N,W] = BSeuCall2Dbasket_RBFFDadap_gs(Nx,g,n,M,Kmul)
%% 2D EU Call RBF-FD with BDF2

funname = 'BSeuCall2Dbasket_RBFFDadap_gs';
datafilename = [funname,'___Nx',num2str(Nx),'_g',num2str(g),'_n',num2str(n),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

names = dir(['./Data/',datafilename]);

if ~isempty(names)
    cd('./Data')
    load(datafilename);
    cd('..')
else
    load('UrefEU.mat')

    tic
    %% Parameters
    phi='gs';
    fit = 3;
    switch n
        case 9
            gamma=0.011; gamma=fit*gamma;
        case 13
            gamma=0.0145; gamma=fit*gamma;
        case 25
            gamma=0.0450; gamma=fit*gamma;
    end

    %% Model
    r=0.03;
    sig1=0.15;
    sig2=0.15;
    rho=0.5;

    T=0.2;
    K=1;

    %% Grid
    Kx=1/Kmul;

    % Nx=100;
    i=1:Nx;
    Ki=2*Kx;
    S=1;

%     g=3; %tune this! 1,2,3,4,5

    c=2*Ki/g;

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

    L=1;
    Nlinsq=sqrt(ceil(N*2-sqrt(N*2)));
    dx=L/(Nlinsq-1);

    % M=100000;
    dt=T/(M-1);
    t=T:-dt:0;

    %% Initial condition
    u0=max(0.5*(xvec+yvec)-Kx,zeros(1,length(xvec)));
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
    % pause()

    %% RBF
    s=[xvec' yvec'];

    % Weights
    indc=findKNearestNeighbors(s,s,n);

    iind=repmat(indin,n,1); iind=iind(:); %n*N
    jind=transpose(indc(indin,:)); jind=jind(:);%n*N
    Wval=zeros(n,numel(indin));  %n*N

    % internal points {
    bb=0;
    lc=zeros(n+1,1);
    for ii=indin
        bb=bb+1;
        %     ii
        %     showsten(1,Nx,xvec,yvec,indc); pause()
        sc=[xvec(ii),yvec(ii)]; xc=sc(:,1); yc=sc(:,2);
        se=s(indc(ii,:),:);

        Rc=xcdist(se,se,1);

        H=Rc(:,:,1);
        hmin=min(min(H(H>0)));
        ep=gamma/hmin;

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
    % } internal points
    Wval=Wval(:);
    W=sparse(iind,jind,Wval,N,N);

    %         display('Weights completed');
    I=speye(size(W));

    %% Integration
    % BDF-1
    u1=u;
    A=I-dt*W;

    b=u1;
    b(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);

    u=A\b;
    u=max(u,zeros(size(u)));

    % BDF-2
    A=I-(2/3)*dt*W;
    rcm=symrcm(A);
    A=A(rcm,rcm);
    [L1, U1]=lu(A);
    for ii=3:M
        u2=u1;
        u1=u;
        b=(4/3)*u1-(1/3)*u2;
        b(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);

        u(rcm)=L1\b(rcm);
        u(rcm)=U1\u(rcm);

        u=max(u,zeros(size(u)));
    end
    tim=toc;
    % figure(3) %solution
    % tri = delaunay(xvec,yvec);
    % trisurf(tri, xvec', yvec', u);
    % axis vis3d
    % axis tight
    % xlabel('S_1');
    % ylabel('S_2');
    % zlabel('V(S_1,S_2)');
    % drawnow;

    %% Error
    indreg=[];
    for ii=1:length(xvec)
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if xvec(ii)>=2/3*Kx && xvec(ii)<=4/3*Kx && yvec(ii)>=2/3*Kx && yvec(ii)<=4/3*Kx
            indreg=[indreg ii];
        end
    end

    xvec=xvec(indreg);
    yvec=yvec(indreg);
    u=u(indreg);

    xvec = K*Kmul*xvec;
    yvec = K*Kmul*yvec;
    dx = K*Kmul*dx;
    u = K*Kmul*u;

    x=[xvec' yvec'];


    uinterp=griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

    err=uinterp'-u;

    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'ep', 'N', 'W', 'Nx', 'g', 'n', 'M', 'Kmul')
    cd('..')
end
