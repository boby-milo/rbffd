function [u,err,tim,x,dx,n,N,W] = BSeuCall2Dbasket_RBFFDadap_smoothed_phs(Nx,g,p,d,nm,M,Kmul)
%% 2D EU Call RBF-FD with BDF2


%% Check if this parameter combination is already computed and if it is just load the result, else compute.
funname = 'BSeuCall2Dbasket_RBFFDadap_smoothed_phs';
datafilename = [funname,'___Nx',num2str(Nx),'_g',num2str(g), '_p',num2str(p),'_d',num2str(d),'_nm',num2str(nm),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

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


    % Initial condition

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

    fu = @(s1, s2) max((1/2)*(s1+s2)-Kx, 0);
    u = smooth4adap(s,dgs,fu,hlocind);

%     u = fu(xvec, yvec);
%     u(indreg) = uind;


    %     fu = @(s1, s2) max((1/2)*(s1+s2)-Kx, 0);
    %     u = fu(xvec, yvec);

    %    indin2=[];
    %    for j=1:Nx
    %        for k=1:Ny
    %            ii=((j-1)*Ny+k);
    %            if (j>1) && (j<Nx) && (k>1) && (k<Ny)
    %                indin2=[indin2 ii];
    %            end
    %        end
    %    end



    %     figure(200)
    %     x = [xvec yvec];
    % tri = delaunay(x(:,1),x(:,2));
    % zz=u-fu(xvec, yvec);
    % figure(99)
    % hold on
    % for i=1:length(x(:,1))
    %     if abs(zz(i))>10^(-14)
    % plot(x(i,1),x(i,2),'.')
    %     end
    % end
    % figure(98)
    % hold on
    % for i=1:length(x(:,1))
    % %    if abs(zz(i))>0
    % plot(x(i,1),x(i,2),'.')
    % %    end
    % end
    % figure(100)
    % trisurf(tri, x(:,1), x(:,2), u-fu(xvec, yvec));
    % shading interp
    % colorbar

    % figure(101)
    % trisurf(tri, x(:,1), x(:,2), u);
    % shading interp
    % colorbar
    % %pause
    % %view(2)
    %
    % u0=max((1/2)*(xvec+yvec)-Kx,zeros(1,length(xvec)));
    % u=u0';

    % figure(1)
    % clf
    % plot(xvec,yvec,'.')
    % hold on
    % plot(xvec(indreg),yvec(indreg),'sq');
    % plot(xvec(indff),yvec(indff),'*')
    % plot(xvec(indcf),yvec(indcf),'^')
    % plot(xvec(indin),yvec(indin),'o')
    % axis equal
    % axis tight
    % hold off
    %
    % figure(2)
    % tri = delaunay(xvec',yvec');
    % trisurf(tri, xvec', yvec', u);
    % shading interp
    % colorbar
    % view(2)
    % axis vis3d
    % % axis equal
    % axis tight
    %
    %     figure(3)
    %      tri = delaunay(xvec',yvec');
    %     trisurf(tri, xvec', yvec', u-fu(xvec,yvec));
    %     shading interp
    %     colorbar
    %     view(2)
    %     axis vis3d
    %     % axis equal
    %     axis tight

    %   pause




    %% Integration
    %     I = speye(size(W));
    %     % BDF-1
    %     u1 = u;
    %     A = I-dt*W;
    %
    %     b = u1;
    %     b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);
    %
    %     u = A\b;
    %     u = max(u,zeros(size(u)));
    %
    %     % BDF-2
    %     A = I-(2/3)*dt*W;
    %     rcm = symrcm(A);
    %     A = A(rcm,rcm);
    %     [L1, U1] = lu(A);
    %     for ii = 3:M
    %         u2 = u1;
    %         u1 = u;
    %         b = (4/3)*u1-(1/3)*u2;
    %         b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);
    %
    %         u(rcm) = L1\b(rcm);
    %         u(rcm) = U1\u(rcm);
    %
    %         u = max(u,zeros(size(u)));
    %     end

    % GMRES
    %  figure(345)
    %   Wsave=W;
    %     save Wsave
    %     pause
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
    b=u1-Abig(indin,indff)*(0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(dt)));

    %     Abig(indin,indff)
    %     pause
    %     xvec(indff)
    %     pause
    %     yvec(indff)
    %     pause
    u = A\b;


    %            usave2II=u2(indin2);
    %    save usave2II
    %   pause

    %    u = max(u,zeros(size(u)));

    % BDF-2
    Abig = I-(2/3)*dt*W;
    A=Abig(indin,indin);
    setup.type='nofill';
    [LA, UA]=ilu(A,setup);
    for ii = 3:M
        t = t+dt;
        u2 = u1;
        u1 = u;



        b = (4/3)*u1-(1/3)*u2 - Abig(indin,indff)*(0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*((ii-1)*dt)));
        %         b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(0.1+(ii-1)*dt));


        [u,~,~,~] = gmres(A,b,200,1e-12,200,LA,UA);





        %      u = max(u,zeros(size(u)));
    end



    ufin(indcf) = 0;
    ufin(indin) = u;
    ufin(indff) = 0.5*(xvec(indff)+yvec(indff)) - Kx*exp(-r*T);
    u = ufin';
    tim = toc;



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
    dx = K*Kmul*dx;
    u = K*Kmul*u(indreg);

    x = [xvec yvec];


    uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

    err = uinterp-u;

    %%
    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'n', 'N', 'W', 'Nx','g', 'p', 'd', 'nm', 'M', 'Kmul')
    cd('..')
    fprintf('...done!\n');
end
end
