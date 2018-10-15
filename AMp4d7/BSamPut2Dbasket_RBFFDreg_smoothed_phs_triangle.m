function [u,err,tim,x,dx,n,N,W,myh,mymax] = BSamPut2Dbasket_RBFFDreg_smoothed_phs_triangle(Nx,p,d,nm,M,Kmul)
%% 2D EU Call RBF-FD with BDF2

funname = 'BSamPut2Dbasket_RBFFDreg_smoothed_phs_triangle';
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
    
    T=1;
    K=1;
    
    %% Grid
    Kx=1/Kmul;
    x=transpose(linspace(0,1,Nx));
    %     dx=x(2)-x(1)
    y=x;
    %
    xvec=[]; yvec=[];
    for ii=1:Nx
        xl=linspace(0,x(ii),ii);
        yl=linspace(x(ii),0,ii);
        %
        xvec=[xvec;xl'];
        yvec=[yvec;yl'];
    end
    
    N=numel(xvec);
    
    ind=1:N;
    indcf=1;
    indff=[length(xvec)-Nx+1:length(xvec)];
    indin=ind; indin([indff,indcf])=[];
    
    L=1;
    Nlinsq=sqrt(ceil(N*2-sqrt(N*2)));
    dx=L/(Nlinsq-1);
    
    
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
    %    n=25;
    n = round(nm*m);
    s = [xvec yvec];
    
    parallel = 0;
    
    %    x
    %    y
    %    pause
    
    [W, ~] = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,x,y,N,n,m,p,indin,phi,d,'reg',parallel,Kx);
    
    fprintf('...done!\n\n')
    % Initial condition
    fprintf('Smoothing...')
    fu = @(s1, s2) max(Kx-(1/2)*(s1+s2), 0);
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
    %    hlocind = hloc(indreg);
    hlocind=dx;
    % %     hlocind
    % %     pause
    xvecind = xvec(indreg);
    yvecind = yvec(indreg);
    uind = smooth4([xvecind, yvecind],fu,hlocind,2);
    
    u = fu(xvec, yvec);
    u(indreg) = uind;
    
    
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
    fprintf(['...done in ', num2str(tim), 's!\n\n'])
    
    fprintf('Saving solution...');
    %% Error
    % indreg = [];
    % for ii = 1:length(xvec)
    %     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
    %     if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
    %         indreg = [indreg ii];
    %     end
    % end
    
    %     xvec=[0;xvec;1];
    %     yvec=[0;yvec;1]
    xvec = K*Kmul*xvec;
    yvec = K*Kmul*yvec;
    %    pause
    dx = K*Kmul*dx;
    u = K*Kmul*u;
    
    x = [xvec yvec];
    
    
    uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');
    %     xulti
    %     pause
    %     yulti
    %     pause
    
    %                 for j=1:Nx-1
    %        for k=1:Ny-1
    %            jj=((j-1)*(Ny-1)+k);
    %            if (j>1) && (j<Nx) && (k>1) && (k<Ny)
    %                indin2=[indin2 jj];
    %            end
    %        end
    %    end
    %    indin2=indin2-1;
    %        usaveII=uinterp(indin2);
    
    %                 for j=1:Nx
    %        for k=1:Ny
    %            jj=((j-1)*(Ny)+k);
    %            if (j>1) && (j<Nx) && (k>1) && (k<Ny)
    %                indin2=[indin2 jj];
    %            end
    %        end
    %    end
    %   % indin2=indin2-1;
    %        usaveII=uinterp(indin2);
    
    %usaveII=uinterp;
    %      usave2II=u2(indin2);
    
    
    err = uinterp-u;
    
    
    %% I am not sure that this numbering can be applied to triangle so I commented this part out
    %     indin2=[];
    %     for j=1:Nx
    %         for k=1:Ny
    %             jj=((j-1)*(Ny)+k);
    %             if (j>1) && (j<Nx) && (k>1) && (k<Ny)
    %                 indin2=[indin2 jj];
    %             end
    %         end
    %     end
    %     %  indin2=indin2-1;
    %     usaveII=err(indin2);
    %         save usaveII
    %   pause
    %%
    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'n', 'N', 'W', 'Nx', 'p', 'd', 'nm', 'M', 'Kmul')
    cd('..')
    fprintf('...done!\n');
end
end
% figure(1)
% tri = delaunay(x(:,1),x(:,2));
% trisurf(tri, x(:,1), x(:,2), (abs(err)));
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
%
% figure(2)
% tri = delaunay(x(:,1),x(:,2));
% trisurf(tri, x(:,1), x(:,2), u);
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
% %caxis([0 0.002])
%
% mymax=0;
% for i=1:length(x(:,1))
%     if x(i,1)>0.8&&x(i,1)<1.2
%         if x(i,2)>0.8&&x(i,2)<1.2
%             if abs(err(i))>mymax
%                 mymax=abs(err(i));
%             end
%         end
%     end
% end
% %myh=1/sqrt(2*length(x(:,1)))
% myh=4/(Nx-1)
% mymax
%
% figure(5)
% tri = delaunay(x(:,1),x(:,2));
% trisurf(tri, x(:,1), x(:,2), (abs(err)));
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
% caxis([0 0.002])
%
% mymax=0;
% for i=1:length(x(:,1))
%     if x(i,1)>0.8&&x(i,1)<1.2
%         if x(i,2)>0.8&&x(i,2)<1.2
%             if abs(err(i))>mymax
%                 mymax=abs(err(i));
%             end
%         end
%     end
% end
% myh=1/sqrt(2*length(x(:,1)))
% mymax
%
%
% figure(15)
% plot(x(:,1),x(:,2),'.')
% axis equal
% axis tight
% end
