function [u,err,tim,x,dx,n,N,W] = BSeuCall2Dbasket_RBFFDadap_smoothed_diag_phs(Nx,g,p,d,nm,M,Kmul)
%% 2D EU Call RBF-FD with BDF2

funname = 'BSeuCall2Dbasket_RBFFDadap_smoothed_diag_phs';
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
    
    T=1;
    K=1;
    
    %% Smoothed Diag
    
    Kmul = 4;
    Kx = 1/Kmul;
    
    i=1:Nx;
    Ki=Kx;
    S=Kmul*Ki;
    
    c=2*Ki/g;
    
    dxi=(1/Nx)*(asinh((S-Ki)/c)-asinh(-Ki/c));
    xi=asinh(-Ki/c)+i*dxi;
    x=[0, Ki+c*sinh(xi)]';
    
    Nx=numel(x);
    
    L=1;
    dx=L/(Nx-1);
    
    % 1D
    fu = @(x) max(x-Kx,0);
    indreg=[];
    for jj = 1:length(x)
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if x(jj) >= 9/10*Kx && x(jj) <= 11/10*Kx
            indreg=[indreg jj];
        end
    end
    
    xind = x(indreg);
    uind = smooth4(xind,fu);
    
    u = fu(x);
    u(indreg)=uind;
    
    % 1D to 2D
    
    dy = diff(x);
    
    xadd = x;
    yadd = x;
    
    xsub = x;
    ysub = x;
    
    xvec = xadd;
    yvec = yadd;
    uvec = u;
    
    for ii = 2:(Nx+1)/2
        xadd = xadd + dy(floor(Nx/2)+ii-1);
        yadd = yadd - dy(floor(Nx/2)+ii-1);
        
        xsub = xsub - dy(floor(Nx/2)+ii-1);
        ysub = ysub + dy(floor(Nx/2)+ii-1);
        
        xaug = [xadd; xsub];
        yaug = [yadd; ysub];
        
        xvec(Nx + (ii-2)*2*Nx + 1 : Nx + (ii-1)*2*Nx) = xaug;
        yvec(Nx + (ii-2)*2*Nx + 1 : Nx + (ii-1)*2*Nx) = yaug;
        
        uvec(Nx + (ii-2)*2*Nx + 1 : Nx + (ii-1)*2*Nx) = [u; u;];
        
    end
    
    br = 1;
    for ii = 1:numel(xvec)
        if xvec(ii) >= -0.000001 && yvec(ii) >= -0.000001 && xvec(ii)+yvec(ii) <= Kx*Kmul+sqrt(2)*dx
            xrot(br) = xvec(ii);
            yrot(br) = yvec(ii);
            urot(br) = uvec(ii);
            
            br = br+1;
            
        end
    end
    
    xvec = xrot;
    yvec = yrot;
    u = urot';
    
    zvec = xvec+yvec;
    maxz = max(zvec);
    
    ind = 1:numel(xvec);
    indff = find(zvec >= maxz - dx);
    indcf = 1;
    
    ind([indff,indcf]) = []; indin = ind;
    
    N = numel(xvec);
    L=1;
    Nlinsq=sqrt(ceil(N*2-sqrt(N*2)));
    dx=L/(Nlinsq-1);
    
    % M=100000;
    dt=T/(M-1);
    t=T:-dt:0;
    
    figure(1)
    clf
    plot(xvec,yvec,'.')
    hold on
    %     plot(xvec(indreg),yvec(indreg),'sq');
    plot(xvec(indff),yvec(indff),'sq')
    plot(xvec(indcf),yvec(indcf),'^')
    plot(xvec(indin),yvec(indin),'o')
    axis equal
    axis tight
    hold off
    
    figure(2)
    tri = delaunay(xvec',yvec');
    trisurf(tri, xvec', yvec', u);
    shading interp
    colorbar
    view(2)
    axis vis3d
    % axis equal
    axis tight
    
    
    %% RBF
    phi='phs';
    
    dim = 2; %problem dimension
    
    m = nchoosek(p+dim, p); %number of polynomial terms;
    n = round(nm*m);
    s = [xvec' yvec'];
    
    parallel = 0;
    [W, ~] = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,N,n,m,p,indin,phi,d,'reg',parallel);
    
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
    % indreg = [];
    % for ii = 1:length(xvec)
    %     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
    %     if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
    %         indreg = [indreg ii];
    %     end
    % end
    
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

figure(10)
tri = delaunay(x(:,1),x(:,2));
trisurf(tri, x(:,1), x(:,2), u);
shading interp
colorbar
view(2)
axis vis3d
% axis equal
axis tight

end
