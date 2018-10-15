function [u,err,tim,x,dx,n,N,W] = BSamPut2Dbasket_RBFFDreg_phs(Nx,p,d,nm,M,Kmul)
%% 2D EU Call RBF-FD with BDF2

funname = 'BSamPut2Dbasket_RBFFDreg_phs';
datafilename = [funname,'___Nx',num2str(Nx),'_p',num2str(p),'_d',num2str(d),'_nm',num2str(nm),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

names = dir(['./Data/',datafilename]);

if ~isempty(names)
    cd('./Data')
    load(datafilename);
    cd('..')
else
    load('UrefAM.mat')

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
    x=transpose(linspace(0,1,Nx));
    % dx=x(2)-x(1)
    y=x;

    xvec=[]; yvec=[];
    for ii=1:Nx
        xl=linspace(0,x(ii),ii);
        yl=linspace(x(ii),0,ii);

        xvec=[xvec,xl];
        yvec=[yvec,yl];
    end

    N=numel(xvec);

    ind=1:N;
    indcf=1;
    indff=[length(xvec)-Nx+1:length(xvec)];
    indin=ind; indin([indff,indcf])=[];

    dx = 1/sqrt(N);

    % M=100000;
    dt=T/(M-1);
    t=T:-dt:0;

    %% Initial condition
    u0=max(Kx-(1/2)*(xvec+yvec),zeros(1,length(xvec)));
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
    phi='phs';

    dim = 2; %problem dimension

    m = nchoosek(p+dim, p); %number of polynomial terms;
    n = round(nm*m);
    s = [xvec' yvec'];

    parallel = 0;
    W = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,x,y,N,n,m,p,indin,phi,d,'reg',parallel,Kx);

   %% Integration
    lambda=zeros(N,1);
    I = speye(size(W));
    % BDF-1
    u1 = u;
    Abig = I-dt*W;
    %    figure(121)
    %     ind
    %     indin
    %     indff
    %     spy(Abig)
    %     pause
    A=Abig(indin,indin);
    
    u1=u1(indin);
    lambda = lambda(indin);
    util=A\(u1+dt*lambda);
    lambdaold=lambda;
    lambda=zeros(numel(indin),1);
    u=util+dt*(lambda-lambdaold);
    for ii=1:numel(indin)
        if u(ii)-(Kx-0.5*(xvec(indin(ii))+yvec(indin(ii))))<0
            u(ii)=Kx-0.5*(xvec(indin(ii))+yvec(indin(ii)));
            lambda(ii)=lambdaold(ii)+(u(ii)-util(ii))/dt;
        end
    end
    
    u=max(u,zeros(size(u)));
    
  
    
    % BDF-2
    Abig=I-(2/3)*dt*W;
    A=Abig(indin,indin);
    setup.type='nofill';
    [LA, UA]=ilu(A,setup);
    
    for ii=3:M
        %     waitbar(ii/M)
        u2=u1;
        u1=u;
        
        b=(4/3)*u1 - (1/3)*u2 + (2/3)*dt*lambda;
        
        %     util(rcm)=L1\b(rcm);
        %     util(rcm)=U1\util(rcm);
        
        [util,~,~,~] = gmres(A,b,200,1e-8,200,LA,UA);
        
        lambdaold=lambda;
        lambda=zeros(numel(indin),1);
        
        u=util+(2/3)*dt*(lambda-lambdaold);
        
        for jj=1:numel(indin)
            if u(jj)-(Kx-0.5*(xvec(indin(jj))+yvec(indin(jj))))<0
                u(jj)=Kx-0.5*(xvec(indin(jj))+yvec(indin(jj)));
                lambda(jj)=lambdaold(jj)+(3/(2*dt))*(u(jj)-util(jj));
            end
        end
        
        u=max(u,zeros(size(u)));
    end
    
    ufin(indcf) = u0(indcf);
    ufin(indin) = u;
    ufin(indff) = zeros(size(indff));
    u = ufin';
    tim = toc;
    
    
%     %% Plotting
%     xvec = K*Kmul*xvec
%     yvec = K*Kmul*yvec
%     dx = K*Kmul*dx
%     u = K*Kmul*u
%     uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');
%     err = uinterp-u;
%     
%     [X Y] = meshgrid(linspace(min(xvec),max(xvec),300), linspace(min(yvec),max(yvec),300));
%     Z = griddata(xvec,yvec,abs(err),X,Y);
%     
%     figure(1)
%     surf(X,Y,Z)
%     axis equal
%     axis tight
%     view(2)
%     drawnow
    
    %% Error
    indreg = [];
    for ii = 1:length(xvec)
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if xvec(ii)>=2/3*Kx && xvec(ii)<=4/3*Kx && yvec(ii)>=2/3*Kx && yvec(ii)<=4/3*Kx
            indreg = [indreg ii];
        end
    end
    
    xvec = K*Kmul*xvec(indreg)';
    yvec = K*Kmul*yvec(indreg)';
    %    pause
    dx = K*Kmul*dx;
    u = K*Kmul*u(indreg);
    
    x = [xvec yvec];
    
    
    uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');
    
    
    
    err = uinterp-u;
    
    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'n', 'N', 'W', 'Nx', 'p', 'd', 'nm', 'M', 'Kmul')
    cd('..')
end
end
