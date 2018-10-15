function [u,err,tim,x,dx,n,N,W] = BSamPut2Dbasket_RBFFDadap_smoothed_phs(Nx,g,p,d,nm,M,Kmul)
%% 2D EU Call RBF-FD with BDF2


%% Check if this parameter combination is already computed and if it is just load the result, else compute.
funname = 'BSamPut2Dbasket_RBFFDadap_smoothed_phs';
datafilename = [funname,'___Nx',num2str(Nx),'_g',num2str(g), '_p',num2str(p),'_d',num2str(d),'_nm',num2str(nm),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

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

    % Nx=100;
    i=1:Nx;
    Ki=2*Kx;
    S=1;

    % g=5; %tune this! 1,2,3,4,5

    c=2*Ki/g;

    dxi=(1/Nx)*(asinh((S-Ki)/c)-asinh(-Ki/c));
    xi=asinh(-Ki/c)+i*dxi;
    x=[0, Ki+c*sinh(xi)];
    y=zeros(numel(x),1);
    Kind=-x(x<=Ki)+Ki;
    
    
    xvec=[]; yvec=[];
    dgind = dsearchn(x',2*Kx);
    dginds = (dgind-3): (dgind+3);
    jj=0;
    for ii=1:numel(x)
        xl=linspace(0,x(ii),ii);
        yl=linspace(x(ii),0,ii);
        
        Nii = numel(xvec);
        if ismember(ii, dginds)
           jj = jj+1;
           indreg = (Nii+1):(Nii+ii);
           dgs{jj} = indreg';
        end
        xvec=[xvec;xl'];
        yvec=[yvec;yl'];
    end

    N=numel(xvec);
    ind=1:N;

    indcf=1;
    indff=(N-numel(x)+1):N;
    indin=ind; indin([indff,indcf])=[];

    dx = 1/sqrt(N);


    dt=T/(M-1);
    t=T:-dt:0;

    %     figure(1)
    %     clf
    %     plot(xvec,yvec,'.')
    %     hold on
    %     plot(xvec(indff),yvec(indff),'sq')
    %     plot(xvec(indcf),yvec(indcf),'^')
    %     plot(xvec(indin),yvec(indin),'ko')
    %     axis equal
    %     axis tight
    %     hold off
    %     pause()


    %% RBF
    phi='phs';

    dim = 2; %problem dimension

    m = nchoosek(p+dim, p); %number of polynomial terms;
    n = round(nm*m);
    s = [xvec yvec];

    parallel = 0;
    [W,hloc] = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,x,y,N,n,m,p,indin,phi,d,'reg',parallel,Kx);


    %% Initial condition

%     fu = @(s1, s2) max((1/2)*(s1+s2)-Kx, 0);
    %
    indreg = [];
    for ii = 1:length(xvec)
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if abs(xvec(ii)+yvec(ii)-2*Kx) < 6*dx
            indreg = [indreg ii];
        end
    end
    %     indreg
    %     pause
    hlocind = hloc(indreg);
    % hlocind=dx;
    % %     hlocind
    % %     pause
%     xvecind = xvec(indreg);
%     yvecind = yvec(indreg);
%     uind = smooth4([xvecind, yvecind],fu,hlocind,2);

    fu = @(s1, s2) max(Kx-(1/2)*(s1+s2), 0);
    u0 = smooth4adap(s,dgs,fu,hlocind);
    u=u0;

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
    
    xvec = K*Kmul*xvec(indreg);
    yvec = K*Kmul*yvec(indreg);
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
