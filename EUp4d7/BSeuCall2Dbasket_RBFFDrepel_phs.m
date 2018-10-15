function [u,err,tim,x,dx,n,N,W] = BSeuCall2Dbasket_RBFFDrepel_phs(Nx,p,d,nm,M,Kmul)
% if nargin == 1
%     pl = 0;
% end
%% Check if this parameter combination is already computed and if it is just load the result, else compute.
funname = 'BSeuCall2Dbasket_RBFFDrepel_phs';
datafilename = [funname,'___Nx',num2str(Nx), '_p',num2str(p),'_d',num2str(d),'_nm',num2str(nm),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

names = dir(['./Data/',datafilename]);

if ~isempty(names)
    cd('./Data')
    load(datafilename);
    cd('..')
else
    load('UrefEU.mat')
    pl = 0;

    %% Parameters
%     M = 100;

    phi = 'phs';
    % ep = 5;
    % p = 4;
    % nm = 5;

%     hlocmul = 1;

    dim = 2; %problem dimension
    m = nchoosek(p+dim, p); %number of polynomial terms;
    n = round(nm*m); n = max(9,n);

    %% Model
    r = 0.03;
    sig1 = 0.15;
    sig2 = 0.15;
    rho = 0.5; %rho=[1 0.5; 0.5 1];

    T = 0.2;
    K = 1;

    %% Grid
    % Kmul = 8;
    Kx = 1/Kmul;

    x_max = Kmul*K;
    y_max = x_max;

    slim = [0, 1, 0, 1];
    xscale = x_max/slim(2);
    yscale = y_max/slim(4);

    [s,N,indin,indcf,indff] = BSeuCall2D_grid_repel(2*Nx, Kx, slim, pl);
    xvec = s(:,1); yvec = s(:,2);
    % Nx = sqrt(N)
    dx = 1/sqrt(N);

%     tri = delaunay(xvec, yvec);
%     trindin = delaunay(xvec(indin), yvec(indin));
%
    tic
    %% Nearest neighbour
    indx = findKNearestNeighbors(s,s,n);

    par = [-r*ones(size(xvec)) ...
        r*xvec ...
        r*yvec ...
        0.5*sig1^2*xvec.^2 ...
        rho*sig1*sig2*xvec.*yvec ...
        0.5*sig2^2*yvec.^2];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [W,wcond,hloc] = RBFFDweights(par,s,N,phi,d,p,n,indin,indx,dim);
%     hloc = hlocmul*hloc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = speye(size(W));
    [wcondmax, icondmax] = max(wcond);
    indcondmax = find(wcond>=1e16);

    %% Integration
    u0 = max(0.5*(xvec+yvec)-Kx, 0);
    u = u0;


    [u,xvec,yvec,Acond] = BSeuCall2D_integrate(u,I,T,M,W,xscale,yscale,xvec,yvec,indin,indff,Kx,r,'gmres');

    tim = toc;

    %% Error
    indreg = [];
    for ii = 1:length(xvec)
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if xvec(ii)>=2/3*K && xvec(ii)<=4/3*K && yvec(ii)>=2/3*K && yvec(ii)<=4/3*K
            indreg = [indreg ii];
        end
    end
    dx = xscale*dx;
    x = [xvec yvec];


    uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

    err = uinterp-u;

    %%
    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'n', 'N', 'W', 'Nx', 'p', 'd', 'nm', 'M', 'Kmul')
    cd('..')

end
end
