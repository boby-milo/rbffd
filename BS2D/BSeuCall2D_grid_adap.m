function [s,N,indres,indin,indcf,indff] = BSeuCall2D_grid_adap(Nx, Kx, sres, slim, pl)

x_min = slim(1);
x_max = slim(2);
y_min = slim(3);
y_max = slim(4);

xres = sres(:,1); yres = sres(:,2);

dx = 1/Nx;

%% Adapted
Kx = 2*Kx;
c = 1/4;
i = 0:Nx;
dxi = (1/Nx) * (asinh((x_max-Kx)/c)-asinh(-Kx/c));
xi = asinh(-Kx/c) + i*dxi;
x = Kx + c*sinh(xi);
x(1) = x_min;

x_max = x(end);

xvec=[]; yvec=[];
for ii = 1:numel(x)
    xl=linspace(0,x(ii),ii)';
    yl=linspace(x(ii),0,ii)';

    xvec=[xvec;xl];
    yvec=[yvec;yl];
end

ind = find(xvec<=x_max & yvec<=y_max & xvec>=x_min & yvec>=y_min);
xvec = xvec(ind);
yvec = yvec(ind);

s = unique([xvec, yvec], 'rows');
xvec = s(:,1);
yvec = s(:,2);

N = numel(xvec);
Nx_res = sqrt(N);

ind = 1:N;
indcf = find(xvec==0 & yvec ==0)';
indff = find(abs(yvec + xvec - x_max) <= 1e-12)'; %, Ny:Ny:N-Ny];
indin = ind;
indin([indff,indcf]) = [];

indres = [];

% indres = [find(xvec == xres(1) & yvec == yres(1)),...
%     find(xvec == xres(2) & yvec == yres(2)),...
%     find(xvec == xres(3) & yvec == yres(3))];

if pl
    figure(1)
    clf
    plot(xvec,yvec,'k.')
    axis equal
    axis tight
    hold on
%     plot(xvec(indin),yvec(indin),'ko')
    plot(xvec(indcf),yvec(indcf),'^')
    plot(xvec(indff),yvec(indff),'sq')
    plot(xvec(indres),yvec(indres),'*')
    drawnow
end
end
