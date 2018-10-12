function [s,N,indres,indin,indcf,indff,inddy0] = HSTeuCall_grid_adap(Nx, Kx, sres, slim, pl)

x_min = slim(1);
x_max = slim(2);
y_min = slim(3);
y_max = slim(4);

% xres = sres(:,1); yres = sres(:,2);

dx = 1/Nx;

%% Adapted
c = 1/8;
Kmul = 1/Kx;
i = 0:Nx;
dxi = (1/Nx) * (asinh((x_max-Kx)/c)-asinh(-Kx/c));
xi = asinh(-Kx/c) + i*dxi;
x = Kx + c*sinh(xi);
x(1) = x_min;

Ny = round(Nx*((y_max-y_min)/x_max));
y = transpose(linspace(y_min,y_max,Ny));

% x = [x; xres];
% y = [y; yres];

[X,Y] = meshgrid(x,y);

xvec = X(:);
yvec = Y(:);

s = unique([xvec, yvec], 'rows');
xvec = s(:,1);
yvec = s(:,2);

N = numel(xvec);
Nx_res = sqrt(N);

ind = 1:N;
indcf = find(xvec==0)';
indff = find(abs(xvec-x_max)<=1e-15)'; %, Ny:Ny:N-Ny];
indin = ind;
indup = find(yvec==y_max)';
inddo = find(yvec==y_min)';

indyup = indup(2:end-1);
indydo = inddo(2:end-1);
inddy0 = [];
% inddy0 = [indyup, indydo];
% inddy0 = indyup;

indin([indff,indcf,inddy0]) = [];
% indin([indff,indcf]) = [];

indres = [];

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
    drawnow
end
end
